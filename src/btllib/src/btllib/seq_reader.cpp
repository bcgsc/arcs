#include "btllib/seq_reader.hpp"
#include "btllib/cstring.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq.hpp"
#include "btllib/status.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <stack>
#include <string>
#include <thread>
#include <vector>

namespace btllib {

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
thread_local std::unique_ptr<decltype(SeqReader::output_queue)::Block>
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  SeqReader::ready_records_array[MAX_SIMULTANEOUS_SEQREADERS];

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
thread_local long
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  SeqReader::ready_records_owners[MAX_SIMULTANEOUS_SEQREADERS] = { 0 };

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
thread_local size_t
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  SeqReader::ready_records_current[MAX_SIMULTANEOUS_SEQREADERS] = { 0 };

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::atomic<long> SeqReader::last_id(0);

SeqReader::SeqReader(const std::string& source_path,
                     const unsigned flags,
                     const unsigned threads)
  : source_path(source_path)
  , source(source_path)
  , flags(flags)
  , threads(threads)
  , buffer_size(short_mode() ? SHORT_MODE_BUFFER_SIZE : LONG_MODE_BUFFER_SIZE)
  , block_size(short_mode() ? SHORT_MODE_BLOCK_SIZE : LONG_MODE_BLOCK_SIZE)
  , cstring_queue(buffer_size, block_size)
  , output_queue(buffer_size, block_size)
  , id(++last_id)
{
  // Parameter sanity check
  check_error(!short_mode() && !long_mode(),
              "SeqReader: no mode selected, either short or long mode flag "
              "must be provided.");
  check_error(short_mode() && long_mode(),
              "SeqReader: short and long mode are mutually exclusive.");
  check_error(threads == 0, "SeqReader: Number of helper threads cannot be 0.");
  start_processors();
  {
    std::unique_lock<std::mutex> lock(format_mutex);
    start_reader();
    format_cv.wait(lock);
  }
}

SeqReader::~SeqReader()
{
  close();
}

void
SeqReader::close() noexcept
{
  bool closed_expected = false;
  if (closed.compare_exchange_strong(closed_expected, true)) {
    try {
      reader_end = true;
      output_queue.close();
      for (auto& pt : processor_threads) {
        pt->join();
      }
      cstring_queue.close();
      reader_thread->join();
      source.close();
    } catch (const std::system_error& e) {
      log_error("SeqReader thread join failure: " + std::string(e.what()));
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }
}

bool
SeqReader::load_buffer()
{
  buffer.start = 0;
  char last = buffer.end > 0 ? buffer.data[buffer.end - 1] : char(0);
  buffer.end = 0;
  do {
    buffer.end += fread(buffer.data.data() + buffer.end,
                        1,
                        buffer.data.size() - buffer.end,
                        source);
  } while (buffer.end < buffer.data.size() && !bool(std::feof(source)));

  if (bool(std::feof(source)) && !buffer.eof_newline_inserted) {
    if (buffer.end < buffer.data.size()) {
      if ((buffer.end == 0 && last != '\n') ||
          (buffer.end > 0 && buffer.data[buffer.end - 1] != '\n')) {
        buffer.data[buffer.end++] = '\n';
      }
      buffer.eof_newline_inserted = true;
    } else if (buffer.data[buffer.data.size() - 1] == '\n') {
      buffer.eof_newline_inserted = true;
    }
    return true;
  }
  return bool(buffer.end);
}

bool
SeqReader::readline_buffer_append(CString& s)
{
  char c = char(0);
  for (; buffer.start < buffer.end && (c = buffer.data[buffer.start]) != '\n';
       ++buffer.start) {
    if (s.s_size >= s.s_cap) {
      s.change_cap(s.s_cap * 2);
    }
    s.s[s.s_size++] = c;
  }
  if (s.s_size >= s.s_cap) {
    s.change_cap(s.s_cap * 2);
  }
  s.s[s.s_size] = '\0';
  if (c == '\n') {
    ++buffer.start;
    return true;
  }
  return false;
}

void
SeqReader::readline_file(CString& s, FILE* f)
{
  s.s_size = getline(&(s.s), &(s.s_cap), f);
}

void
SeqReader::readline_file_append(CString& s, FILE* f)
{
  readline_file(tmp, f);
  if (s.s_size + tmp.s_size + 1 > s.s_cap) {
    s.change_cap(s.s_size + tmp.s_size + 1);
  }
  memcpy(s.s + s.s_size, tmp.s, tmp.s_size + 1);
  s.s_size += tmp.s_size;
}

bool
SeqReader::file_at_end(FILE* f)
{
  if (std::ferror(f) == 0) {
    auto c = std::fgetc(f);
    if (c == EOF) {
      return true;
    }
    std::ungetc(c, f);
    return false;
  }
  return true;
}

int
SeqReader::getc_buffer()
{
  if (buffer.start < buffer.end) {
    return buffer.data[buffer.start++];
  }
  return EOF;
}

int
SeqReader::ungetc_buffer(const int c)
{
  if (buffer.start > 0) {
    buffer.data[--buffer.start] = char(c);
    return c;
  }
  return EOF;
}

void
SeqReader::update_cstring_records(OrderQueueSPMC<RecordCString>::Block& records,
                                  size_t& counter)
{
  records.count++;
  if (records.count == block_size) {
    records.num = counter++;
    cstring_queue.write(records);
    records.num = 0;
    records.count = 0;
  }
}

#define BTLLIB_SEQREADER_FORMAT_CHECK(INPUT_FORMAT, READ_MODULE)               \
  if ((++module_counter, (READ_MODULE).buffer_valid(buf, bufsize))) {          \
    format = Format::INPUT_FORMAT;                                             \
    module_in_use = module_counter;                                            \
  } else

#define BTLLIB_SEQREADER_FORMAT_READ(READ_MODULE)                              \
  if (module_in_use == ++module_counter) {                                     \
    read_from_buffer(READ_MODULE, records, counter);                           \
    read_transition(READ_MODULE, records, counter);                            \
    read_from_file(READ_MODULE, records, counter);                             \
  } else

void
SeqReader::determine_format()
{
  load_buffer();
  bool empty = buffer.end - buffer.start == 1;
  check_warning(empty, std::string(source_path) + " is empty.");

  if (empty) {
    return;
  }

  auto* const buf = buffer.data.data() + buffer.start;
  const auto bufsize = buffer.end - buffer.start;

  int module_counter = 0;

  BTLLIB_SEQREADER_FORMAT_CHECK(FASTA, fasta_module)
  BTLLIB_SEQREADER_FORMAT_CHECK(FASTA, multiline_fasta_module)
  BTLLIB_SEQREADER_FORMAT_CHECK(FASTQ, fastq_module)
  BTLLIB_SEQREADER_FORMAT_CHECK(FASTQ, multiline_fastq_module)
  BTLLIB_SEQREADER_FORMAT_CHECK(SAM, sam_module)
  BTLLIB_SEQREADER_FORMAT_CHECK(GFA2, gfa2_module)
  {
    format = Format::INVALID;
    log_error(std::string(source_path) + " source file is in invalid format!");
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

void
SeqReader::start_reader()
{
  reader_thread = std::unique_ptr<std::thread>(new std::thread([this]() {
    {
      std::unique_lock<std::mutex> lock(format_mutex);
      determine_format();
      format_cv.notify_all();
    }

    size_t counter = 0;
    decltype(cstring_queue)::Block records(block_size);

    if (get_format() != Format::UNDETERMINED) {
      int module_counter = 0;

      BTLLIB_SEQREADER_FORMAT_READ(fasta_module)
      BTLLIB_SEQREADER_FORMAT_READ(multiline_fasta_module)
      BTLLIB_SEQREADER_FORMAT_READ(fastq_module)
      BTLLIB_SEQREADER_FORMAT_READ(multiline_fastq_module)
      BTLLIB_SEQREADER_FORMAT_READ(sam_module)
      BTLLIB_SEQREADER_FORMAT_READ(gfa2_module)
      {
        log_error("SeqReader: No reading module was enabled.");
        std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    }

    reader_end = true;
    if (records.count > 0) {
      records.num = counter++;
      cstring_queue.write(records);
    }
    for (unsigned i = 0; i < threads; i++) {
      if (i == 0) {
        dummy_block_num = counter;
      }
      decltype(cstring_queue)::Block dummy(block_size);
      dummy.num = counter++;
      dummy.count = 0;
      cstring_queue.write(dummy);
    }
  }));
}

#undef BTLLIB_SEQREADER_FORMAT_CHECK
#undef BTLLIB_SEQREADER_FORMAT_READ

void
SeqReader::start_processors()
{
  processor_threads.reserve(threads);
  for (unsigned i = 0; i < threads; i++) {
    processor_threads.push_back(
      std::unique_ptr<std::thread>(new std::thread([this]() {
        decltype(cstring_queue)::Block records_in(block_size);
        decltype(output_queue)::Block records_out(block_size);
        for (;;) {
          cstring_queue.read(records_in);
          for (size_t i = 0; i < records_in.count; i++) {
            records_out.data[i].seq = std::string(
              records_in.data[i].seq, records_in.data[i].seq.size());
            auto& seq = records_out.data[i].seq;
            rtrim(seq);

            records_out.data[i].qual = std::string(
              records_in.data[i].qual, records_in.data[i].qual.size());
            auto& qual = records_out.data[i].qual;
            rtrim(qual);

            char *first_whitespace = nullptr, *last_whitespace = nullptr;
            for (size_t j = 0; j < records_in.data[i].header.size(); j++) {
              if (bool(std::isspace(records_in.data[i].header[j]))) {
                if (first_whitespace == nullptr) {
                  first_whitespace = records_in.data[i].header + j;
                }
                last_whitespace = records_in.data[i].header + j;
              } else if (last_whitespace != nullptr) {
                break;
              }
            }
            size_t id_start = (format == Format::FASTA ||
                               format == Format::FASTQ || format == Format::SAM)
                                ? 1
                                : 0;

            switch (format) {
              case Format::FASTA:
                check_error(records_in.data[i].header.empty(),
                            "SeqReader: Invalid FASTA header");
                check_error(records_in.data[i].header[0] != '>',
                            "SeqReader: Unexpected character in a FASTA file.");
                break;
              case Format::FASTQ:
                check_error(records_in.data[i].header.empty(),
                            "SeqReader: Invalid FASTQ header");
                check_error(records_in.data[i].header[0] != '@',
                            "SeqReader: Unexpected character in a FASTQ file.");
                break;
              default:
                break;
            }

            if (first_whitespace == nullptr) {
              records_out.data[i].id =
                std::string(records_in.data[i].header + id_start,
                            records_in.data[i].header.size() - id_start);
              records_out.data[i].comment = "";
            } else {
              records_out.data[i].id = std::string(
                records_in.data[i].header + id_start,
                first_whitespace - records_in.data[i].header - id_start);
              records_out.data[i].comment = std::string(
                last_whitespace + 1,
                records_in.data[i].header.size() -
                  (last_whitespace - records_in.data[i].header) - 1);
            }
            records_in.data[i].header.clear();

            auto& id = records_out.data[i].id;
            auto& comment = records_out.data[i].comment;
            rtrim(id);
            rtrim(comment);

            if (trim_masked()) {
              const auto len = seq.length();
              size_t trim_start = 0, trim_end = seq.length();
              while (trim_start <= len && bool(islower(seq[trim_start]))) {
                trim_start++;
              }
              while (trim_end > 0 && bool(islower(seq[trim_end - 1]))) {
                trim_end--;
              }
              seq.erase(trim_end);
              seq.erase(0, trim_start);
              if (!qual.empty()) {
                qual.erase(trim_end);
                qual.erase(0, trim_start);
              }
            }
            if (fold_case()) {
              for (auto& c : seq) {
                char old = c;
                c = CAPITALS[(unsigned char)(c)];
                if (!bool(c)) {
                  log_error(std::string("A sequence contains invalid "
                                        "IUPAC character: ") +
                            old);
                  std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
                }
              }
            }
            records_out.data[i].num = records_in.num * block_size + i;
          }
          records_out.count = records_in.count;
          records_out.num = records_in.num;
          if (records_out.count == 0) {
            if (records_out.num == dummy_block_num) {
              output_queue.write(records_out);
            }
            break;
          }
          output_queue.write(records_out);
        }
      })));
  }
}

SeqReader::Record
SeqReader::read()
{
  if (ready_records_owners[id % MAX_SIMULTANEOUS_SEQREADERS] != id) {
    ready_records_array[id % MAX_SIMULTANEOUS_SEQREADERS] =
      std::unique_ptr<decltype(output_queue)::Block>(
        new decltype(output_queue)::Block(block_size));
    ready_records_owners[id % MAX_SIMULTANEOUS_SEQREADERS] = id;
    ready_records_current[id % MAX_SIMULTANEOUS_SEQREADERS] = 0;
  }
  auto& ready_records =
    *(ready_records_array[id % MAX_SIMULTANEOUS_SEQREADERS]);
  auto& current = ready_records_current[id % MAX_SIMULTANEOUS_SEQREADERS];
  if (current >= ready_records.count) {
    ready_records.count = 0;
    output_queue.read(ready_records);
    if (ready_records.count == 0) {
      close();
      ready_records = decltype(output_queue)::Block(block_size);
      return Record();
    }
    current = 0;
  }
  return std::move(ready_records.data[current++]);
}

OrderQueueMPMC<SeqReader::Record>::Block
SeqReader::read_block()
{
  decltype(SeqReader::read_block()) block(block_size);
  output_queue.read(block);
  if (block.count == 0) {
    close();
  }
  return block;
}

} // namespace btllib
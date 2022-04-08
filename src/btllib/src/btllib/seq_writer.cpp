#include "btllib/seq_writer.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/seq.hpp"
#include "btllib/status.hpp"

#include <cstdio>
#include <mutex>
#include <string>

namespace btllib {

SeqWriter::SeqWriter(const std::string& sink_path, Format format, bool append)
  : sink_path(sink_path)
  , sink(sink_path, append)
  , closed(false)
  , format(format)
  , headerchar(format == FASTA ? '>' : '@')
{}

void
SeqWriter::close()
{
  if (!closed) {
    sink.close();
    closed = true;
  }
}

void
SeqWriter::write(const std::string& id,
                 const std::string& comment,
                 const std::string& seq,
                 const std::string& qual)
{
  check_error(seq.empty(), "Attempted to write empty sequence.");
  for (const auto& c : seq) {
    if (!bool(COMPLEMENTS[(unsigned char)(c)])) {
      log_error(std::string("A sequence contains invalid IUPAC character: ") +
                c);
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }

  std::string output;
  output.reserve(1 + id.size() + 1 + comment.size() + 1 + seq.size() + 3 +
                 qual.size() + 1);
  output += headerchar;
  if (!id.empty()) {
    output += id;
  }
  if (!comment.empty()) {
    output += " ";
    output += comment;
  }
  output += '\n';

  output += seq;
  output += '\n';

  if (format == FASTQ) {
    check_error(seq.size() != qual.size(),
                "Quality must be the same length as sequence.");
    output += "+\n";
    output += qual;
    output += '\n';
  }

  {
    std::unique_lock<std::mutex> lock(mutex);
    if (fwrite(output.c_str(), 1, output.size(), sink) != output.size()) {
      log_error("SeqWriter: fwrite failed: " + get_strerror());
    }
  }
}

} // namespace btllib
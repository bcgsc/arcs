#include "btllib/indexlr.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/status.hpp"

#include "btllib_config.hpp"

#include <cassert>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

const static std::string PROGNAME = "indexlr";
const static std::string VERSION = btllib::PROJECT_VERSION;
const static size_t OUTPUT_PERIOD_SHORT = 512;
const static size_t OUTPUT_PERIOD_LONG = 2;
const static size_t INITIAL_OUTPUT_STREAM_SIZE = 100;
const static size_t QUEUE_SIZE = 64;
const static size_t MAX_THREADS = 5;
const static size_t DEFAULT_THREADS = MAX_THREADS;

static void
print_error_msg(const std::string& msg)
{
  std::cerr << PROGNAME << ' ' << VERSION << ": " << msg << std::endl;
}

static void
print_usage()
{
  std::cerr << "Usage: " << PROGNAME
            << " -k K -w W [-r repeat_bf_path] [-s solid_bf_path] [--id] "
               "[--bx] [--pos] [--seq] "
               "[-o FILE] FILE...\n\n"
               "  -k K        Use K as k-mer size.\n"
               "  -w W        Use W as sliding-window size.\n"
               "  --id        Include input sequence ids in the output. "
               "(Default if --bx is not "
               "provided)\n"
               "  --bx        Include input sequence barcodes in the output.\n"
               "  --len       Include input sequence length in the output.\n"
               "  --pos       Include minimizer positions in the output "
               "(appended with : after "
               "minimizer value).\n"
               "  --strand    Include minimizer strands in the output "
               "(appended with : after minimizer "
               "value).\n"
               "  --seq       Include minimizer sequences in the output "
               "(appended with : after "
               "minimizer value).\n"
               "              If a combination of --pos, --strand, and --seq "
               "options are provided, "
               "they're appended in the --pos, --strand, --seq order after the "
               "minimizer value.\n"
               "  --long      Enable long mode which is more efficient for "
               "long sequences (e.g. long "
               "reads, contigs, reference).\n"
               "  -r repeat_bf_path  Use a Bloom filter to filter out "
               "repetitive minimizers.\n"
               "  -s solid_bf_path  Use a Bloom filter to only select solid "
               "minimizers.\n"
               "  -o FILE     Write output to FILE, default is stdout.\n"
               "  -t T        Use T number of threads (default 5, max 5) per "
               "input file.\n"
               "  -v          Show verbose output.\n"
               "  --help      Display this help and exit.\n"
               "  --version   Display version and exit.\n"
               "  FILE        Space separated list of FASTA/Q files."
            << std::endl;
}

int
main(int argc, char* argv[])
{
  int c;
  int optindex = 0;
  int help = 0, version = 0;
  bool verbose = false;
  unsigned k = 0, w = 0, t = DEFAULT_THREADS;
  bool w_set = false;
  bool k_set = false;
  int with_id = 0, with_bx = 0, with_len = 0, with_pos = 0, with_strand = 0,
      with_seq = 0;
  std::unique_ptr<btllib::KmerBloomFilter> repeat_bf, solid_bf;
  bool with_repeat = false, with_solid = false;
  int long_mode = 0;
  std::string outfile("-");
  bool failed = false;
  static const struct option longopts[] = {
    { "id", no_argument, &with_id, 1 },
    { "bx", no_argument, &with_bx, 1 },
    { "len", no_argument, &with_len, 1 },
    { "pos", no_argument, &with_pos, 1 },
    { "strand", no_argument, &with_strand, 1 },
    { "seq", no_argument, &with_seq, 1 },
    { "long", no_argument, &long_mode, 1 },
    { "help", no_argument, &help, 1 },
    { "version", no_argument, &version, 1 },
    { nullptr, 0, nullptr, 0 }
  };
  while ((c = getopt_long(argc, // NOLINT(concurrency-mt-unsafe)
                          argv,
                          "k:w:o:t:vr:s:",
                          longopts,
                          &optindex)) != -1) {
    switch (c) {
      case 0:
        break;
      case 'k':
        k_set = true;
        k = std::stoul(optarg);
        break;
      case 'w':
        w_set = true;
        w = std::stoul(optarg);
        break;
      case 'o':
        outfile = optarg;
        break;
      case 't':
        t = std::stoul(optarg);
        break;
      case 'v':
        verbose = true;
        break;
      case 'r': {
        with_repeat = true;
        std::cerr << "Loading repeat Bloom filter from " << optarg << std::endl;
        try {
          repeat_bf = std::unique_ptr<btllib::KmerBloomFilter>(
            new btllib::KmerBloomFilter(optarg));
        } catch (const std::exception& e) {
          std::cerr << e.what() << '\n';
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
        }
        std::cerr << "Finished loading repeat Bloom filter" << std::endl;
        break;
      }
      case 's': {
        with_solid = true;
        std::cerr << "Loading solid Bloom filter from " << optarg << std::endl;
        try {
          solid_bf = std::unique_ptr<btllib::KmerBloomFilter>(
            new btllib::KmerBloomFilter(optarg));
        } catch (const std::exception& e) {
          std::cerr << e.what() << '\n';
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
        }
        std::cerr << "Finished loading solid Bloom filter" << std::endl;
        break;
      }
      default:
        std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }
  if (t > MAX_THREADS) {
    t = MAX_THREADS;
    std::cerr << (PROGNAME + ' ' + VERSION + ": Using more than " +
                  std::to_string(MAX_THREADS) +
                  " threads does not scale, reverting to 5.\n")
              << std::flush;
  }
  std::vector<std::string> infiles(&argv[optind], &argv[argc]);
  if (argc < 2) {
    print_usage();
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  if (help != 0) {
    print_usage();
    std::exit(EXIT_SUCCESS); // NOLINT(concurrency-mt-unsafe)
  } else if (version != 0) {
    std::cerr << PROGNAME << ' ' << VERSION << std::endl;
    std::exit(EXIT_SUCCESS); // NOLINT(concurrency-mt-unsafe)
  }
  if (!k_set) {
    print_error_msg("missing option -- 'k'");
    failed = true;
  } else if (k == 0) {
    print_error_msg("option has incorrect value -- 'k'");
    failed = true;
  }
  if (!w_set) {
    print_error_msg("missing option -- 'w'");
    failed = true;
  } else if (w == 0) {
    print_error_msg("option has incorrect value -- 'w'");
    failed = true;
  }
  if (infiles.empty()) {
    print_error_msg("missing file operand");
    failed = true;
  }
  if (failed) {
    std::cerr << "Try '" << PROGNAME << " --help' for more information.\n";
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }

  unsigned flags = 0;
  if (bool(with_bx) && !bool(with_id)) {
    flags |= btllib::Indexlr::Flag::NO_ID;
  }
  if (bool(with_bx)) {
    flags |= btllib::Indexlr::Flag::BX;
  }
  if (bool(with_seq)) {
    flags |= btllib::Indexlr::Flag::SEQ;
  }
  if (bool(long_mode)) {
    flags |= btllib::Indexlr::Flag::LONG_MODE;
  } else {
    flags |= btllib::Indexlr::Flag::SHORT_MODE;
  }

  btllib::Indexlr::Record record;
  FILE* out;
  if (outfile == "-") {
    out = stdout;
  } else {
#ifdef __linux__
    out = fopen(outfile.c_str(), "we");
#else
    out = fopen(outfile.c_str(), "w"); // NOLINT(android-cloexec-fopen)
#endif
  }
  for (auto& infile : infiles) {
    std::unique_ptr<btllib::Indexlr> indexlr;
    try {
      if (with_repeat && with_solid) {
        flags |= btllib::Indexlr::Flag::FILTER_IN;
        flags |= btllib::Indexlr::Flag::FILTER_OUT;
        indexlr = std::unique_ptr<btllib::Indexlr>(
          new btllib::Indexlr(infile,
                              k,
                              w,
                              flags,
                              t,
                              verbose,
                              solid_bf->get_bloom_filter(),
                              repeat_bf->get_bloom_filter()));
      } else if (with_repeat) {
        flags |= btllib::Indexlr::Flag::FILTER_OUT;
        indexlr = std::unique_ptr<btllib::Indexlr>(new btllib::Indexlr(
          infile, k, w, flags, t, verbose, repeat_bf->get_bloom_filter()));
      } else if (with_solid) {
        flags |= btllib::Indexlr::Flag::FILTER_IN;
        indexlr = std::unique_ptr<btllib::Indexlr>(new btllib::Indexlr(
          infile, k, w, flags, t, verbose, solid_bf->get_bloom_filter()));
      } else {
        indexlr = std::unique_ptr<btllib::Indexlr>(
          new btllib::Indexlr(infile, k, w, flags, t, verbose));
      }
    } catch (const std::exception& e) {
      std::cerr << e.what() << '\n';
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
    std::queue<std::string> output_queue;
    std::mutex output_queue_mutex;
    std::condition_variable queue_empty, queue_full;
    size_t max_seen_output_size = INITIAL_OUTPUT_STREAM_SIZE;
    const size_t output_period =
      bool(long_mode) ? OUTPUT_PERIOD_LONG : OUTPUT_PERIOD_SHORT;
    try {
      std::unique_ptr<std::thread> info_compiler(new std::thread([&]() {
        std::stringstream ss;
        while ((record = indexlr->read())) {
          if (bool(with_id) || (!bool(with_id) && !bool(with_bx))) {
            ss << record.id << '\t';
          }
          if (bool(with_bx)) {
            ss << record.barcode << '\t';
          }
          if (bool(with_len)) {
            ss << record.readlen << '\t';
          }
          int j = 0;
          for (const auto& min : record.minimizers) {
            if (j > 0) {
              ss << ' ';
            }
            ss << min.out_hash;
            if (bool(with_pos)) {
              ss << ':' << min.pos;
            }
            if (bool(with_strand)) {
              ss << ':' << (min.forward ? '+' : '-');
            }
            if (bool(with_seq)) {
              ss << ':' << min.seq;
            }
            j++;
          }
          ss << '\n';
          if (record.num % output_period == output_period - 1) {
            auto ss_str = ss.str();
            max_seen_output_size =
              std::max(max_seen_output_size, ss_str.size());
            std::unique_lock<std::mutex> lock(output_queue_mutex);
            while (output_queue.size() == QUEUE_SIZE) {
              queue_full.wait(lock);
            }
            output_queue.push(std::move(ss_str));
            queue_empty.notify_one();
            lock.unlock();
            std::string newstring;
            newstring.reserve(max_seen_output_size);
            ss.str(newstring);
          }
        }
        {
          std::unique_lock<std::mutex> lock(output_queue_mutex);
          if (!ss.str().empty()) {
            output_queue.push(ss.str());
          }
          output_queue.push(std::string());
          queue_empty.notify_one();
        }
      }));
      std::unique_ptr<std::thread> output_worker(new std::thread([&]() {
        std::string to_write;
        for (;;) {
          {
            std::unique_lock<std::mutex> lock(output_queue_mutex);
            while (output_queue.empty()) {
              queue_empty.wait(lock);
            }
            to_write = std::move(output_queue.front());
            output_queue.pop();
            queue_full.notify_one();
          }
          if (to_write.empty()) {
            break;
          }
          btllib::check_error(
            fwrite(to_write.c_str(), 1, to_write.size(), out) !=
              to_write.size(),
            "Indexlr: fwrite failed.");
        }
      }));
      info_compiler->join();
      output_worker->join();
    } catch (const std::exception& e) {
      std::cerr << e.what() << '\n';
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }
  if (out != stdout) {
    fclose(out);
  }

  return 0;
}

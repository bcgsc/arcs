#include "btllib/bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/status.hpp"

#include "cpptoml.h"

#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <sys/stat.h>

namespace btllib {

static unsigned
pop_cnt_byte(uint8_t x)
{
  return ((0x876543210 >>                                              // NOLINT
           (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >> // NOLINT
          ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) &  // NOLINT
         0xf;                                                          // NOLINT
}

BloomFilter::BloomFilter(size_t bytes, unsigned hash_num, std::string hash_fn)
  : bytes(
      size_t(std::ceil(double(bytes) / sizeof(uint64_t)) * sizeof(uint64_t)))
  , array_size(get_bytes() / sizeof(array[0]))
  , array_bits(array_size * CHAR_BIT)
  , hash_num(hash_num)
  , hash_fn(std::move(hash_fn))
  , array(new std::atomic<uint8_t>[array_size])
{
  // Parameter sanity check
  check_error(bytes == 0, "BloomFilter: memory budget must be >0!");
  check_error(hash_num == 0, "BloomFilter: number of hash values must be >0!");
  check_error(hash_num > MAX_HASH_VALUES,
              "BloomFilter: number of hash values cannot be over 1024!");
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  std::memset((void*)array.get(), 0, array_size * sizeof(array[0]));
}

void
BloomFilter::insert(const uint64_t* hashes)
{
  for (unsigned i = 0; i < hash_num; ++i) {
    const auto normalized = hashes[i] % array_bits;
    array[normalized / CHAR_BIT] |= BIT_MASKS[normalized % CHAR_BIT];
  }
}

bool
BloomFilter::contains(const uint64_t* hashes) const
{
  for (unsigned i = 0; i < hash_num; ++i) {
    const auto normalized = hashes[i] % array_bits;
    const auto mask = BIT_MASKS[normalized % CHAR_BIT];
    if (!bool(array[normalized / CHAR_BIT] & mask)) {
      return false;
    }
  }
  return true;
}

bool
BloomFilter::contains_insert(const uint64_t* hashes)
{
  uint8_t found = 1;
  for (unsigned i = 0; i < hash_num; ++i) {
    const auto normalized = hashes[i] % array_bits;
    const auto bitpos = normalized % CHAR_BIT;
    const auto mask = BIT_MASKS[bitpos];
    found &= ((array[normalized / CHAR_BIT].fetch_or(mask) >> bitpos) & 1);
  }
  return bool(found);
}

uint64_t
BloomFilter::get_pop_cnt() const
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for default(none) reduction(+ : pop_cnt)
  for (size_t i = 0; i < array_size; ++i) {
    pop_cnt += pop_cnt_byte(array[i]);
  }
  return pop_cnt;
}

double
BloomFilter::get_occupancy() const
{
  return double(get_pop_cnt()) / double(array_bits);
}

double
BloomFilter::get_fpr() const
{
  return std::pow(get_occupancy(), double(hash_num));
}

bool
BloomFilterInitializer::check_file_signature(
  std::ifstream& ifs,
  const std::string& expected_signature,
  std::string& file_signature)
{
  std::getline(ifs, file_signature);
  return file_signature == expected_signature;
}

std::shared_ptr<cpptoml::table>
BloomFilterInitializer::parse_header(const std::string& expected_signature)
{
  struct stat buffer
  {};
  btllib::check_error(stat(path.c_str(), &buffer) != 0,
                      "BloomFilterInitializer: " + get_strerror() + ": " +
                        path);
  btllib::check_error(ifs.fail(),
                      "BloomFilterInitializer: failed to open " + path);

  std::string file_signature;
  if (!check_file_signature(ifs, expected_signature, file_signature)) {
    log_error(std::string("File signature does not match (possibly version "
                          "mismatch) for file:\n") +
              path + '\n' + "Expected signature:\t" + expected_signature +
              '\n' + "File signature:    \t" + file_signature);
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }

  /* Read bloom filter line by line until it sees "[HeaderEnd]"
  which is used to mark the end of the header section and
  assigns the header to a char array*/
  std::string toml_buffer(file_signature + '\n');
  std::string line;
  bool header_end_found = false;
  while (bool(std::getline(ifs, line))) {
    toml_buffer.append(line + '\n');
    if (line == "[HeaderEnd]") {
      header_end_found = true;
      break;
    }
  }
  if (!header_end_found) {
    log_error("Pre-built bloom filter does not have the correct header end.");
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES; i++) {
    std::getline(ifs, line);
  }

  // Send the char array to a stringstream for the cpptoml parser to parse
  std::istringstream toml_stream(toml_buffer);
  cpptoml::parser toml_parser(toml_stream);
  const auto header_config = toml_parser.parse();

  // Obtain header values from toml parser and assign them to class members
  const auto header_string =
    file_signature.substr(1, file_signature.size() - 2); // Remove [ ]
  return header_config->get_table(header_string);
}

BloomFilter::BloomFilter(const std::string& path)
  : BloomFilter::BloomFilter(
      std::make_shared<BloomFilterInitializer>(path, BLOOM_FILTER_SIGNATURE))
{}

BloomFilter::BloomFilter(const std::shared_ptr<BloomFilterInitializer>& bfi)
  : bytes(*(bfi->table->get_as<decltype(bytes)>("bytes")))
  , array_size(bytes / sizeof(array[0]))
  , array_bits(array_size * CHAR_BIT)
  , hash_num(*(bfi->table->get_as<decltype(hash_num)>("hash_num")))
  , hash_fn(bfi->table->contains("hash_fn")
              ? *(bfi->table->get_as<decltype(hash_fn)>("hash_fn"))
              : "")
  , array(new std::atomic<uint8_t>[array_size])
{
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  bfi->ifs.read((char*)array.get(),
                std::streamsize(array_size * sizeof(array[0])));
}

void
BloomFilter::save(const std::string& path,
                  const cpptoml::table& table,
                  const char* data,
                  const size_t n)
{
  std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);

  ofs << table << "[HeaderEnd]\n";
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES; i++) {
    if (i == 1) {
      ofs << "  <binary data>";
    }
    ofs << '\n';
  }

  ofs.write(data, std::streamsize(n));
}

bool
BloomFilter::check_file_signature(const std::string& path,
                                  const std::string& signature)
{
  std::ifstream ifs(path);
  std::string file_signature;
  return BloomFilterInitializer::check_file_signature(
    ifs, signature, file_signature);
}

void
BloomFilter::save(const std::string& path)
{
  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  if (!hash_fn.empty()) {
    header->insert("hash_fn", get_hash_fn());
  }
  std::string header_string = BLOOM_FILTER_SIGNATURE;
  header_string =
    header_string.substr(1, header_string.size() - 2); // Remove [ ]
  root->insert(header_string, header);
  save(path, *root, (char*)array.get(), array_size * sizeof(array[0]));
}

KmerBloomFilter::KmerBloomFilter(size_t bytes, unsigned hash_num, unsigned k)
  : k(k)
  , bloom_filter(bytes, hash_num, HASH_FN)
{}

void
KmerBloomFilter::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    bloom_filter.insert(nthash.hashes());
  }
}

unsigned
KmerBloomFilter::contains(const char* seq, size_t seq_len) const
{
  unsigned count = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    if (bloom_filter.contains(nthash.hashes())) {
      count++;
    }
  }
  return count;
}

unsigned
KmerBloomFilter::contains_insert(const char* seq, size_t seq_len)
{
  unsigned count = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    if (bloom_filter.contains_insert(nthash.hashes())) {
      count++;
    }
  }
  return count;
}

KmerBloomFilter::KmerBloomFilter(const std::string& path)
  : KmerBloomFilter::KmerBloomFilter(
      std::make_shared<BloomFilterInitializer>(path,
                                               KMER_BLOOM_FILTER_SIGNATURE))
{}

KmerBloomFilter::KmerBloomFilter(
  const std::shared_ptr<BloomFilterInitializer>& bfi)
  : k(*(bfi->table->get_as<decltype(k)>("k")))
  , bloom_filter(bfi)
{
  check_error(bloom_filter.hash_fn != HASH_FN,
              "KmerBloomFilter: loaded hash function (" + bloom_filter.hash_fn +
                ") is different from the one used by default (" + HASH_FN +
                ").");
}

void
KmerBloomFilter::save(const std::string& path)
{
  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  header->insert("hash_fn", get_hash_fn());
  header->insert("k", get_k());
  std::string header_string = KMER_BLOOM_FILTER_SIGNATURE;
  header_string =
    header_string.substr(1, header_string.size() - 2); // Remove [ ]
  root->insert(header_string, header);

  BloomFilter::save(path,
                    *root,
                    (char*)bloom_filter.array.get(),
                    bloom_filter.array_size * sizeof(bloom_filter.array[0]));
}

SeedBloomFilter::SeedBloomFilter(size_t bytes,
                                 unsigned k,
                                 const std::vector<std::string>& seeds,
                                 unsigned hash_num_per_seed)
  : seeds(seeds)
  , parsed_seeds(parse_seeds(seeds))
  , kmer_bloom_filter(bytes, hash_num_per_seed, k)
{
  for (const auto& seed : seeds) {
    check_error(k != seed.size(),
                "SeedBloomFilter: passed k (" + std::to_string(k) +
                  ") not equal to passed spaced seed size (" +
                  std::to_string(seed.size()) + ")");
  }
}

void
SeedBloomFilter::insert(const char* seq, size_t seq_len)
{
  SeedNtHash nthash(
    seq, seq_len, parsed_seeds, get_hash_num_per_seed(), get_k());
  while (nthash.roll()) {
    for (size_t s = 0; s < seeds.size(); s++) {
      kmer_bloom_filter.bloom_filter.insert(nthash.hashes() +
                                            s * get_hash_num_per_seed());
    }
  }
}

std::vector<std::vector<unsigned>>
SeedBloomFilter::contains(const char* seq, size_t seq_len) const
{
  std::vector<std::vector<unsigned>> hit_seeds;
  SeedNtHash nthash(
    seq, seq_len, parsed_seeds, get_hash_num_per_seed(), get_k());
  while (nthash.roll()) {
    hit_seeds.emplace_back();
    for (size_t s = 0; s < seeds.size(); s++) {
      if (kmer_bloom_filter.bloom_filter.contains(
            nthash.hashes() + s * get_hash_num_per_seed())) {
        hit_seeds.back().push_back(s);
      }
    }
  }
  return hit_seeds;
}

std::vector<std::vector<unsigned>>
SeedBloomFilter::contains_insert(const char* seq, size_t seq_len)
{
  std::vector<std::vector<unsigned>> hit_seeds;
  SeedNtHash nthash(
    seq, seq_len, parsed_seeds, get_hash_num_per_seed(), get_k());
  while (nthash.roll()) {
    hit_seeds.emplace_back();
    for (size_t s = 0; s < seeds.size(); s++) {
      if (kmer_bloom_filter.bloom_filter.contains_insert(
            nthash.hashes() + s * get_hash_num_per_seed())) {
        hit_seeds.back().push_back(s);
      }
    }
  }
  return hit_seeds;
}

double
SeedBloomFilter::get_fpr() const
{
  const double single_seed_fpr =
    std::pow(get_occupancy(), get_hash_num_per_seed());
  return 1 - std::pow(1 - single_seed_fpr, seeds.size());
}

SeedBloomFilter::SeedBloomFilter(const std::string& path)
  : SeedBloomFilter::SeedBloomFilter(
      std::make_shared<BloomFilterInitializer>(path,
                                               SEED_BLOOM_FILTER_SIGNATURE))
{}

SeedBloomFilter::SeedBloomFilter(
  const std::shared_ptr<BloomFilterInitializer>& bfi)
  : seeds(*(bfi->table->get_array_of<std::string>("seeds")))
  , parsed_seeds(parse_seeds(seeds))
  , kmer_bloom_filter(bfi)
{}

void
SeedBloomFilter::save(const std::string& path)
{
  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  header->insert("hash_fn", get_hash_fn());
  header->insert("k", get_k());
  auto seeds_array = cpptoml::make_array();
  for (const auto& seed : seeds) {
    seeds_array->push_back(seed);
  }
  header->insert("seeds", seeds_array);
  std::string header_string = SEED_BLOOM_FILTER_SIGNATURE;
  header_string =
    header_string.substr(1, header_string.size() - 2); // Remove [ ]
  root->insert(header_string, header);

  BloomFilter::save(path,
                    *root,
                    (char*)kmer_bloom_filter.bloom_filter.array.get(),
                    kmer_bloom_filter.bloom_filter.array_size *
                      sizeof(kmer_bloom_filter.bloom_filter.array[0]));
}

} // namespace btllib

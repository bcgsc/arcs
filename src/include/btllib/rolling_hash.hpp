#ifndef BTLLIB_ROLLING_HASH_HPP
#define BTLLIB_ROLLING_HASH_HPP

#include "nthash.hpp"

#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace btllib
{

class RollingHash;
class SeedRollingHash;
using SpacedSeed = std::vector<unsigned>;
static std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string> &seed_strings);

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */
class RollingHash
{

public:
  /**
   * Constructor.
   * @param seq DNA sequence to be hashed
   * @param seq_len length of seq
   * @param k k-mer size
   * @param hash_num number of hashes
   */
  RollingHash(const char *seq, size_t seq_len, unsigned k, unsigned hash_num);

  /**
   * Constructor.
   * @param seq DNA sequence to be hashed
   * @param k k-mer size
   * @param hash_num number of hashes
   */
  RollingHash(const std::string &seq, unsigned k, unsigned hash_num);

  /**
   * Calculate the next hash value
   * @return true on success and false otherwise
   */
  bool roll();

  const uint64_t *hashes() const;

  size_t get_pos() const { return pos; }
  unsigned get_k() const { return k; }
  unsigned get_hash_num() const { return hash_num; }

protected:
  /** Initialize internal state of iterator */
  bool init();

  const char *seq;
  const size_t seq_len;
  const unsigned k;
  const unsigned hash_num;
  size_t pos = 0;
  std::vector<uint64_t> hashes_vector;
  uint64_t forward_hash = 0;
  uint64_t reverse_hash = 0;
};

class SeedRollingHash : public RollingHash
{

public:
  SeedRollingHash(const char *seq,
                  size_t seq_len,
                  unsigned k,
                  const std::vector<SpacedSeed> &seeds,
                  unsigned hash_num_per_seed);
  SeedRollingHash(const std::string &seq,
                  unsigned k,
                  const std::vector<SpacedSeed> &seeds,
                  unsigned hash_num_per_seed);
  SeedRollingHash(const char *seq,
                  size_t seq_len,
                  unsigned k,
                  const std::vector<std::string> &seeds,
                  unsigned hash_num_per_seed);
  SeedRollingHash(const std::string &seq,
                  unsigned k,
                  const std::vector<std::string> &seeds,
                  unsigned hash_num_per_seed);

  unsigned get_hash_num_per_seed() const { return hash_num_per_seed; }

  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> hash_components;

  bool roll()
  {
    init();
    hashes_vector[0] = 0;
    return true;
  }

private:
  bool init()
  {

    if (k > seq_len)
    {
      pos = std::numeric_limits<std::size_t>::max();
      return false;
    }
    std::cerr << "checkpoint1" << std::endl;

    unsigned max_block_length = 0;
    for (unsigned i = 0; i < seeds.size(); ++i)
    {
      const btllib::SpacedSeed &seed = seeds.at(i);
      std::vector<std::pair<unsigned, unsigned>> one_blocks;
      unsigned block_length = 0;
      if (seed.at(0) > max_block_length)
      {
        max_block_length = seed.at(0);
      }
      one_blocks.emplace_back(std::make_pair(0, seed.at(0)));

      for (unsigned j = 1; j < seed.size(); ++j)
      {
        block_length = seed.at(j) - seed.at(j - 1) - 1;
        if (seed.at(j) - seed.at(j - 1) > 1)
        {
          one_blocks.emplace_back(std::make_pair(seed.at(j - 1) + 1, block_length));
        }
        if (block_length > max_block_length)
        {
          max_block_length = seed.at(0);
        }
        /*std::cerr << "checkpoint11" << std::endl;

        std::cerr << seed.at(j) << std::endl;
        if (seed.at(j) == 1)
        {
          std::cerr << "checkpoint10" << std::endl;

          ++block_length;
          if (at_ones == false)
          {
            at_ones = true;
            start_of_block = j;
          }
        }
        if (seed.at(j) == 0 && at_ones == true)
        {
          at_ones = false;
          one_blocks.emplace_back(std::make_pair(start_of_block, block_length));
          if (block_length > max_block_length)
          {
            max_block_length = block_length;
          }
          block_length = 0;
        }*/
      }
      block_length = seq_len - seed.back() - 1;
      one_blocks.emplace_back(std::make_pair(seed.back() + 1, block_length));

      //one_blocks.emplace_back(std::make_pair(start_of_block, block_length));
      if (block_length > max_block_length)
      {
        max_block_length = block_length;
      }
      one_blocks_of_seeds.emplace_back(one_blocks);
    }
    std::cerr << max_block_length << std::endl;
    std::cerr << "checkpoint2" << std::endl;

    hash_components.resize(max_block_length + 1);
    for (const auto &one_blocks_of_seed : one_blocks_of_seeds)
    {
      for (const auto &one_block : one_blocks_of_seed)
      {
        auto &kmer_size = std::get<1>(one_block);
        std::cerr << kmer_size << std::endl;
        if (hash_components[kmer_size].size() == 0)
        {

          pos = 0;
          hash_components[kmer_size] = std::vector<std::pair<uint64_t, uint64_t>>();
          unsigned posN = 0;
          while ((pos < seq_len - kmer_size + 1) && !(NTC64(seq + pos, kmer_size, posN, hash_components[kmer_size])))
          {
            pos += posN + 1;
          }
          if (pos > seq_len - kmer_size)
          {
            pos = std::numeric_limits<std::size_t>::max();
            return false;
          }
          ++pos;

          while (pos < seq_len - kmer_size + 1)
          {
            if (seed_tab[(unsigned char)(seq[pos + kmer_size - 1])] == seedN)
            {
              pos += kmer_size;
              for (unsigned i = 0; i < kmer_size; ++i)
              {
                hash_components[kmer_size].emplace_back(std::make_pair((uint64_t)0, (uint64_t)0));
              }
              while ((pos < seq_len - kmer_size + 1) && !(NTC64(seq + pos, kmer_size, posN, hash_components[kmer_size])))
              {
                pos += posN + 1;
                hash_components[kmer_size].emplace_back(std::make_pair((uint64_t)0, (uint64_t)0));
              }
            }
            NTMC64(seq[pos - 1], seq[pos - 1 + kmer_size], kmer_size, hash_components[kmer_size]);
            ++pos;
          }
        }
      }
    }
    return true;
  }

  const unsigned hash_num_per_seed;
  std::vector<SpacedSeed> seeds;
  std::vector<std::vector<std::pair<unsigned, unsigned>>> one_blocks_of_seeds;
};

inline RollingHash::RollingHash(const char *seq,
                                size_t seq_len,
                                unsigned k,
                                unsigned hash_num)
    : seq(seq), seq_len(seq_len), k(k), hash_num(hash_num)
{
  hashes_vector.resize(hash_num);
}

inline RollingHash::RollingHash(const std::string &seq,
                                unsigned k,
                                unsigned hash_num)
    : RollingHash(seq.c_str(), seq.size(), k, hash_num)
{
}

inline SeedRollingHash::SeedRollingHash(const char *seq,
                                        size_t seq_len,
                                        unsigned k,
                                        const std::vector<SpacedSeed> &seeds,
                                        unsigned hash_num_per_seed)
    : RollingHash(seq, seq_len, k, seeds.size() * hash_num_per_seed), hash_num_per_seed(hash_num_per_seed), seeds(seeds)
{
}

inline SeedRollingHash::SeedRollingHash(const std::string &seq,
                                        unsigned k,
                                        const std::vector<SpacedSeed> &seeds,
                                        unsigned hash_num_per_seed)
    : RollingHash(seq, k, seeds.size() * hash_num_per_seed), hash_num_per_seed(hash_num_per_seed), seeds(seeds)
{
}

inline SeedRollingHash::SeedRollingHash(const char *seq,
                                        size_t seq_len,
                                        unsigned k,
                                        const std::vector<std::string> &seeds,
                                        unsigned hash_num_per_seed)
    : RollingHash(seq, seq_len, k, seeds.size() * hash_num_per_seed), hash_num_per_seed(hash_num_per_seed), seeds(parse_seeds(seeds))
{
}

inline SeedRollingHash::SeedRollingHash(const std::string &seq,
                                        unsigned k,
                                        const std::vector<std::string> &seeds,
                                        unsigned hash_num_per_seed)
    : RollingHash(seq, k, seeds.size() * hash_num_per_seed), hash_num_per_seed(hash_num_per_seed), seeds(parse_seeds(seeds))
{
}

static std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string> &seed_strings)
{
  std::vector<SpacedSeed> seed_set;
  for (const auto &seed_string : seed_strings)
  {
    SpacedSeed seed;
    size_t pos = 0;
    for (const auto &c : seed_string)
    {
      if (c != '1')
      {
        seed.push_back(pos);
      }
      ++pos;
    }
    seed_set.push_back(seed);
  }
  return seed_set;
}

// NOLINTNEXTLINE
#define ROLLING_HASH_INIT(CLASS, NTHASH_CALL)         \
  inline bool CLASS::init()                           \
  {                                                   \
    if (k > seq_len)                                  \
    {                                                 \
      pos = std::numeric_limits<std::size_t>::max();  \
      return false;                                   \
    }                                                 \
    unsigned posN = 0;                                \
    while ((pos < seq_len - k + 1) && !(NTHASH_CALL)) \
    {                                                 \
      pos += posN + 1;                                \
    }                                                 \
    if (pos > seq_len - k)                            \
    {                                                 \
      pos = std::numeric_limits<std::size_t>::max();  \
      return false;                                   \
    }                                                 \
    ++pos;                                            \
    return true;                                      \
  }

// NOLINTNEXTLINE
#define ROLLING_HASH_ROLL(CLASS, NTHASH_CALL)                 \
  inline bool CLASS::roll()                                   \
  {                                                           \
    if (pos == 0)                                             \
    {                                                         \
      return init();                                          \
    }                                                         \
    if (pos > seq_len - k)                                    \
    {                                                         \
      return false;                                           \
    }                                                         \
    if (seed_tab[(unsigned char)(seq[pos + k - 1])] == seedN) \
    {                                                         \
      pos += k;                                               \
      return init();                                          \
    }                                                         \
    (NTHASH_CALL);                                            \
    ++pos;                                                    \
    return true;                                              \
  }

ROLLING_HASH_INIT(RollingHash,
                  NTMC64(seq + pos,
                         k,
                         hash_num,
                         forward_hash,
                         reverse_hash,
                         posN,
                         hashes_vector.data()))
ROLLING_HASH_ROLL(RollingHash,
                  NTMC64(seq[pos - 1],
                         seq[pos - 1 + k],
                         k,
                         hash_num,
                         forward_hash,
                         reverse_hash,
                         hashes_vector.data()))

#undef ROLLING_HASH_INIT
#undef ROLLING_HASH_ROLL

inline const uint64_t *
RollingHash::hashes() const
{
  return hashes_vector.data();
}

} // namespace btllib

#endif

#include "btllib/counting_bloom_filter.hpp"

#include "helpers.hpp"

#include <cstdio>
#include <iostream>
#include <string>

int
main()
{
  std::cerr << "Testing CountingBloomFilter" << std::endl;
  btllib::CountingBloomFilter8 cbf(1024 * 1024, 3);

  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 100, 200, 300 });

  TEST_ASSERT_EQ(cbf.contains({ 1, 10, 100 }), 2);
  TEST_ASSERT_EQ(cbf.contains({ 100, 200, 300 }), 1);
  TEST_ASSERT_EQ(cbf.contains({ 1, 20, 100 }), 0);

  auto filename = get_random_name(64);
  cbf.save(filename);

  btllib::CountingBloomFilter8::is_bloom_file(filename);
  btllib::CountingBloomFilter8 cbf2(filename);

  TEST_ASSERT_EQ(cbf2.contains({ 1, 10, 100 }), 2);
  TEST_ASSERT_EQ(cbf2.contains({ 100, 200, 300 }), 1);
  TEST_ASSERT_EQ(cbf2.contains({ 1, 20, 100 }), 0);

  TEST_ASSERT_EQ(cbf2.contains_insert({ 9, 99, 999 }), 0);
  TEST_ASSERT_EQ(cbf2.contains_insert({ 9, 99, 999 }), 1);

  TEST_ASSERT_EQ(cbf2.insert_contains({ 9, 99, 999 }), 3);
  TEST_ASSERT_EQ(cbf2.insert_contains({ 9, 99, 999 }), 4);

  TEST_ASSERT_EQ(cbf2.contains_insert_thresh({ 9, 99, 999 }, 5), 4);
  TEST_ASSERT_EQ(cbf2.contains_insert_thresh({ 9, 99, 999 }, 5), 5);
  TEST_ASSERT_EQ(cbf2.contains_insert_thresh({ 9, 99, 999 }, 5), 5);

  TEST_ASSERT_EQ(cbf2.insert_thresh_contains({ 9, 99, 999 }, 6), 6);
  TEST_ASSERT_EQ(cbf2.insert_thresh_contains({ 9, 99, 999 }, 6), 6);

  std::remove(filename.c_str());

  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  TEST_ASSERT_EQ(seq.size(), seq2.size());

  std::cerr << "Testing KmerCountingBloomFilter" << std::endl;
  btllib::KmerCountingBloomFilter8 kbf(1024 * 1024, 4, seq.size() / 2);
  kbf.insert(seq);
  const auto expected = seq.size() - seq.size() / 2 + 1;
  TEST_ASSERT_EQ(kbf.contains(seq), expected);
  TEST_ASSERT_LE(kbf.contains(seq2), 1);

  filename = get_random_name(64);
  kbf.save(filename);

  btllib::KmerCountingBloomFilter8 kbf2(filename);
  TEST_ASSERT_EQ(kbf2.contains(seq), expected);
  TEST_ASSERT_LE(kbf2.contains(seq2), 1);

  TEST_ASSERT_EQ(kbf2.contains_insert(seq), expected);
  TEST_ASSERT_EQ(kbf2.contains_insert(seq), expected * 2);

  TEST_ASSERT_EQ(kbf2.insert_contains(seq), expected * 4);
  TEST_ASSERT_EQ(kbf2.insert_contains(seq), expected * 5);

  TEST_ASSERT_EQ(kbf2.contains_insert_thresh(seq, 6), expected * 5);
  TEST_ASSERT_EQ(kbf2.contains_insert_thresh(seq, 6), expected * 6);
  TEST_ASSERT_EQ(kbf2.contains_insert_thresh(seq, 6), expected * 6);

  TEST_ASSERT_EQ(kbf2.insert_thresh_contains(seq, 7), expected * 7);
  TEST_ASSERT_EQ(kbf2.insert_thresh_contains(seq, 7), expected * 7);

  std::remove(filename.c_str());

  std::cerr << "Testing KmerCountingBloomfilter with multiple threads"
            << std::endl;

  std::vector<std::string> present_seqs;
  std::vector<unsigned> inserts;
  std::vector<std::string> absent_seqs;
  for (size_t i = 0; i < 100; i++) {
    present_seqs.push_back(get_random_seq(100));
    inserts.push_back(get_random(1, 10));
    absent_seqs.push_back(get_random_seq(100));
  }
  std::vector<std::string> present_seqs2 = present_seqs;

  btllib::KmerCountingBloomFilter8 kbf_multithreads(100 * 1024 * 1024, 4, 100);
#pragma omp parallel shared(                                                   \
  present_seqs, present_seqs2, inserts, absent_seqs, kbf_multithreads)
  {
    while (true) {
      std::string seq;
      unsigned insert_count;
      bool end = false;
#pragma omp critical
      {
        if (present_seqs.empty()) {
          end = true;
        } else {
          seq = present_seqs.back();
          insert_count = inserts.back();
          present_seqs.pop_back();
          inserts.pop_back();
        }
      }
      if (end) {
        break;
      }
      for (size_t i = 0; i < insert_count; i++) {
        kbf_multithreads.insert(seq);
      }
    }
  }

  unsigned false_positives = 0;
#pragma omp parallel shared(present_seqs,                                      \
                            present_seqs2,                                     \
                            inserts,                                           \
                            absent_seqs,                                       \
                            kbf_multithreads)                                  \
                            reduction(+:false_positives)
  {
    while (true) {
      std::string seq;
      bool end = false;
#pragma omp critical
      {
        if (absent_seqs.empty()) {
          end = true;
        } else {
          seq = absent_seqs.back();
          absent_seqs.pop_back();
        }
      }
      if (end) {
        break;
      }
      false_positives += kbf_multithreads.contains(seq);
    }
  }
  std::cerr << "False positives = " << false_positives << std::endl;
  TEST_ASSERT_LT(false_positives, 10);

  int more_than_1 = 0;
#pragma omp parallel shared(present_seqs,                                      \
                            present_seqs2,                                     \
                            inserts,                                           \
                            absent_seqs,                                       \
                            kbf_multithreads)                                  \
                            reduction(+:more_than_1)
  {
    while (true) {
      std::string seq;
      bool end = false;
#pragma omp critical
      {
        if (present_seqs2.empty()) {
          end = true;
        } else {
          seq = present_seqs2.back();
          present_seqs2.pop_back();
        }
      }
      if (end) {
        break;
      }
      if (kbf_multithreads.contains(seq) > 1) {
        more_than_1++;
      }
      TEST_ASSERT(kbf_multithreads.contains(seq));
    }
  }
  std::cerr << "Seqs with more than 1 presence = " << more_than_1 << std::endl;
  TEST_ASSERT_GT(more_than_1, 5);

  return 0;
}
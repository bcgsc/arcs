#include "../include/btllib/seq_reader.hpp"

#include <iostream>

int
main()
{
  for (const auto& record : btllib::SeqReader(
         "my_reads.fq.gz", btllib::SeqReader::Flag::SHORT_MODE)) {
    std::cout << record.seq << '\n' << record.qual << '\n';
  }

  /*
   * In Python, this would be:
   *
   * with btllib.SeqReader("my_reads.fq.gz", btllib.SeqReaderFlag.SHORT_MODE) as
   reader: for record in reader: print(record.seq, record.qual)
   *
   */

  return 0;
}
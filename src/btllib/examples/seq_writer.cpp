#include "../include/btllib/seq_writer.hpp"

#include <iostream>

int
main()
{
  btllib::SeqWriter writer("my_reads.fq.gz", btllib::SeqWriter::FASTQ);
  writer.write("read1", "comment1", "ACGATCAGCAGCATG", "!@$&!#$&@(#*@)%");
  writer.write("read2", "comment2", "GTACG", "&%(@(");
  return 0;
}
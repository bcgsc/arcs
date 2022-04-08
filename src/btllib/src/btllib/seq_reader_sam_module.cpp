#include "btllib/seq_reader_sam_module.hpp"
#include "btllib/cstring.hpp"
#include "btllib/process_pipeline.hpp"
#include "btllib/seq.hpp"
#include "btllib/status.hpp"

#include <cstdio>
#include <memory>
#include <thread>

namespace btllib {

bool
SeqReaderSamModule::buffer_valid(const char* buffer, const size_t size)
{
  enum Column
  {
    QNAME = 1,
    FLAG,
    RNAME,
    POS,
    MAPQ,
    CIGAR,
    RNEXT,
    PNEXT,
    TLEN,
    SEQ,
    QUAL
  };

  size_t current = 0;

  while (current < size && buffer[current] == '@') {
    while (current < size && buffer[current] != '\n') {
      current++;
    }
    current++;
  }

  Column column = QNAME;
  unsigned char c;
  while (current < size) {
    c = buffer[current];
    if (c == '\n') {
      break;
    }
    if (c == '\t') {
      if (current > 0 && !bool(std::isspace(buffer[current - 1]))) {
        column = Column(int(column) + 1);
      } else {
        return false;
      }
    } else {
      switch (column) {
        case QNAME:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case FLAG:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case RNAME:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case POS:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case MAPQ:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case CIGAR:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case RNEXT:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case PNEXT:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case TLEN:
          if (!(bool(std::isdigit(c)) || c == '-')) {
            return false;
          }
          break;
        case SEQ:
          if (!bool(COMPLEMENTS[c])) {
            return false;
          }
          break;
        case QUAL:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        default:
          break;
      }
    }
    current++;
  }

  return current >= size || column >= QUAL;
}

} // namespace btllib

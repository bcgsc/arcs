#include "btllib/seq_reader_multiline_fasta_module.hpp"
#include "btllib/seq.hpp"

namespace btllib {

bool
SeqReaderMultilineFastaModule::buffer_valid(const char* buffer,
                                            const size_t size)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ,
    IN_TRANSITION
  };
  State state = IN_HEADER_1;
  while (current < size) {
    c = buffer[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '>') {
          state = IN_HEADER_2;
        } else {
          return false;
        }
        break;
      case IN_HEADER_2:
        if (c == '\n') {
          state = IN_SEQ;
        }
        break;
      case IN_SEQ:
        if (c == '\n') {
          state = IN_TRANSITION;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
      case IN_TRANSITION:
        if (c == '>') {
          state = IN_HEADER_2;
          break;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        state = IN_SEQ;
        break;
    }
    current++;
  }
  return true;
}

} // namespace btllib
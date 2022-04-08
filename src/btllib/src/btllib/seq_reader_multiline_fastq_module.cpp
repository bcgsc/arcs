#include "btllib/seq_reader_multiline_fastq_module.hpp"
#include "btllib/seq.hpp"

namespace btllib {

bool
SeqReaderMultilineFastqModule::buffer_valid(const char* buffer,
                                            const size_t size)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ,
    IN_TRANSITION,
    IN_PLUS_2,
    IN_QUAL
  };
  size_t seqlen = 0, quallen = 0;
  State state = IN_HEADER_1;
  while (current < size) {
    c = buffer[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '@') {
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
        } else if (c != '\r') {
          if (!bool(COMPLEMENTS[c])) {
            return false;
          }
          seqlen++;
        }
        break;
      case IN_TRANSITION:
        if (c == '+') {
          state = IN_PLUS_2;
          break;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        seqlen++;
        state = IN_SEQ;
        break;
      case IN_PLUS_2:
        if (c == '\n') {
          state = IN_QUAL;
        }
        break;
      case IN_QUAL:
        if (quallen < seqlen) {
          if (c != '\r' && c != '\n') {
            if (c < '!' || c > '~') {
              return false;
            }
            quallen++;
          }
        } else if (c == '\n') {
          seqlen = 0;
          quallen = 0;
          state = IN_HEADER_1;
        } else if (c != '\r') {
          return false;
        }
        break;
    }
    current++;
  }
  return true;
}

} // namespace btllib
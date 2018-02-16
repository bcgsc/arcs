#ifndef _PARSE_UTIL_H_
#define _PARSE_UTIL_H_ 1

/**
 * Extract the first Chromium barcode sequence from a given string
 * (FASTQ comment or list of SAM tags). The barcode is assumed
 * to be embedded with the format BX:Z:<BARCODE_SEQ>-1, as
 * per "longranger" FASTQ/SAM output.
 */
static inline std::string getBarcodeSeq(const std::string& str)
{
    size_t start = str.find("BX:Z:");
    if (start == std::string::npos) {
        return std::string();
    }

    /* find next whitespace or EOL after "BX:Z:" */
    size_t end = str.find_first_of(" \t\r\n\v", start);
    if (end == std::string::npos) {
        end = str.length();
    }

    /* sanity check: word must end with "-1" */
    if (str.substr(end - 2, 2) != "-1") {
        return std::string();
    }

    /* adjust `start` and `end` to barcode sequence boundaries */
    start += 5;
    end -= 2;

    /* sanity check: Chromium barcode seq should be >= 16 bases long */
    if (end - start < 16) {
        return std::string();
    }

    return str.substr(start, end - start);
}

#endif

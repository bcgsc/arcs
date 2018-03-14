#ifndef PARSE_UTIL_H
#define PARSE_UTIL_H 1

#include "Common/StringUtil.h"
#include <string>

/**
 * Extract the first Chromium barcode sequence from a given string
 * (FASTQ comment or list of SAM tags). The barcode is expected
 * to match the format BX:Z:<BARCODE_SEQ>.
 */
static inline std::string parseBarcode(const std::string& str)
{
    size_t start = str.find("BX:Z:");
    if (start == std::string::npos)
        return std::string();
    start += 5;

    // Find the next whitespace or EOL after "BX:Z:".
    size_t end = str.find_first_of(" \t\r\n", start);
    if (end == std::string::npos)
        end = str.length();

    return str.substr(start, end - start);
}

#endif

#ifndef SAM_H
#define SAM_H 1

#include <string>

/** Extract the specified SAM tag from a string.
 * @param tag the SAM tag, including two colons, for example "BX:Z:"
 */
template <size_t N>
static inline std::string parseSAMTag(const std::string& s, const char (&tag)[N])
{
    size_t start = s.find(tag);
    if (start == std::string::npos)
        return std::string();
    start += N - 1;

    // Find the next whitespace or EOL after "BX:Z:".
    size_t end = s.find_first_of(" \t\r\n", start);
    if (end == std::string::npos)
        end = s.length();

    return s.substr(start, end - start);
}

/**
 * Extract the first Chromium barcode sequence from a given string
 * (FASTQ comment or list of SAM tags). The barcode is expected
 * to match the format BX:Z:<BARCODE>.
 */
static inline std::string parseBXTag(const std::string& s)
{
    return parseSAMTag(s, "BX:Z:");
}

#endif

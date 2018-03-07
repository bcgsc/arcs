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

/*
 * Check if SAM flag is one of the accepted ones.
 */
bool checkFlag(int flag) {
    return (flag == 99 || flag == 163 || flag == 83 || flag == 147);
}

/*
 * Check if character is one of the accepted ones.
 */
bool checkChar(char c) {
    return (c == 'M' || c == '=' || c == 'X' || c == 'I');
}

/* Returns true if seqence only contains ATGC and is of length indexLen */
bool checkIndex(std::string seq) {
    for (int i = 0; i < static_cast<int>(seq.length()); i++) {
        char c = toupper(seq[i]);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C')
            return false;
    }
    //return (static_cast<int>(seq.length()) == params.indexLen);
    return true;
}

/*
 * Calculate the sequence identity from the cigar string
 * sequence length, and tags.
 */
double calcSequenceIdentity(const std::string& line, const std::string& cigar, const std::string& seq) {

    int qalen = 0;
    std::stringstream ss;
    for (auto i = cigar.begin(); i != cigar.end(); ++i) {
        if (!isdigit(*i)) {
            if (checkChar(*i)) {
                ss << "\t";
                int value = 0;
                ss >> value;
                qalen += value;
                ss.str("");
            } else {
                ss.str("");
            }
        } else {
            ss << *i;
        }
    }

    int edit_dist = 0;
    std::size_t found = line.find("NM:i:");
    if (found!=std::string::npos) {
        edit_dist = std::strtol(&line[found + 5], 0, 10);
    }

    double si = 0;
    if (qalen != 0) {
        double mins = qalen - edit_dist;
        double div = mins/seq.length();
        si = div * 100;
    }

    return si;
}


#endif

#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H 1

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
	extern bool colourSpace;
	extern int verbose;
	extern unsigned streakThreshold;
	extern unsigned threads;
	extern int rank;
}

#endif

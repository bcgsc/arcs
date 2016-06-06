/*
 * Dynamicofstream.h
 *	A output stream that wraps around common output stream
 *	objects, but with the added benefit of dynamically using
 *	the appropriate type depending on output filename
 *
 *	Uses normal ofstream in all cases except:
 *	Gzip (.gz) file extensions uses zlib (gzstream)
 *
 *  Created on: Jun 19, 2013
 *      Author: cjustin
 */
//@TODO: Thread so that gzipping process does not run on main thread?

#ifndef DYNAMICOFSTREAM_H_
#define DYNAMICOFSTREAM_H_

#include <string>
#include <stdint.h>

using namespace std;

class Dynamicofstream{
public:
	Dynamicofstream(const string &filename);
//	void write(const string &input);
//	Dynamicofstream& operator <<(Dynamicofstream& out, const string& o);
	ostream& operator <<(const string& o);
	ostream& operator <<(unsigned o);
	void close();
	virtual ~Dynamicofstream();
private:
	ostream* filestream;

	//@TODO: Not happy with having to store this like this
	//Should figure out better way and refactor code
	bool gz;

};

#endif /* DYNAMICOFSTREAM_H_ */

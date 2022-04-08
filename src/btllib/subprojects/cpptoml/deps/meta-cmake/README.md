# meta-cmake: A collection of cmake helper scripts
The scripts here are used in the [MeTA toolkit][meta]'s build system.

The most useful ones are probably `FindLIBCXX.cmake` and
`FindOrBuildICU.cmake`.

- `FindLIBCXX.cmake` is a standard CMake module for finding libc++ and a
   suitable ABI (libc++abi or libcxxrt). It also provides a macro
   `set_libcxx_required_flags()` to set CMAKE_REQUIRED_FLAGS for performing
   compilation checks using libc++.

- `FindOrBuildICU.cmake` is used to either locate a suitable ICU library
   (matching a specified version exactly) or build one using
   `ExternalProject_Add`.

- `CompilerKludges.cmake` provides a macro for setting a bunch of
  preprocessor variables that are used as workarounds for compiler bugs,
  standard library bugs, or missing standard library features.

- `SetClangOptions.cmake` provides a macro that is probably only useful in
   the context of MeTA.

# License

Copyright (c) 2016 Chase Geigle

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[meta]: https://meta-toolkit.org

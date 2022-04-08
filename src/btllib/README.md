[Bioinformatics Technology Lab](http://www.birollab.ca/) common code library in C++ with Python wrappers.

[![Build Status](https://dev.azure.com/bcgsc/btl_public/_apis/build/status/bcgsc.btllib)](https://dev.azure.com/bcgsc/btl_public/_build/latest?definitionId=1)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/bcgsc/btllib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/bcgsc/btllib/context:cpp)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/bcgsc/btllib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/bcgsc/btllib/alerts/)

Platforms
---
- Linux
- MacOS

Documentation
---
[Docs page](https://bcgsc.github.io/btllib/)

Download
---
The recommended way is to download the [latest release](https://github.com/bcgsc/btllib/releases/latest).

Dependencies
---
- Build
  * GCC 6+ or Clang 5+ with OpenMP
  * Python 3.5+
  * Meson and Ninja Python3 packages, CMake (If you are a user and not a developer, these will be automatically installed to a temporary directory)
- Run time
  * SAMtools for reading SAM, BAM, and CRAM files.
  * gzip, tar, pigz, bzip2, xz, lrzip, zip, and/or 7zip for compressing/decompressing files. Not all of these are necessary, only the ones whose compressions you'll be using. 
  * wget for downloading sequences from a URL.

For users
---
- Copy the root `btllib` directory into your project
- Run `btllib/compile`
- C++
  * Link your code with `btllib/install/lib/libbtllib.a` (pass `-L /path/to/btllib/install/lib -l btllib` flags to the compiler).
  * `#include` any header from the `btllib/install/include` directory (pass `-I /path/to/btllib/install/include` flag to the compiler).
  * `btllib` uses `C++11` features, so that standard should be enabled at a minimum.
- Python wrappers
  * The wrappers correspond one-to-one with C++ code so any functions and classes can be used under the same name. The only exception are nested classes which are prefixed with outer class name (e.g. `btllib::SeqReader::Flag` in C++ versus `btllib.SeqReaderFlag` in Python).
  * Use `PYTHONPATH` environment variable or `sys.path.append()` in your Python code to include `/path/to/btllib/install/python` directory
  * Include the library with `import btllib`
- Executables
  * btllib generated executables can be found in `/path/to/btllib/install/bin` directory. Append that path to the `PATH` environment variable to make it available to your shell.

For developers
---
- Initial setup:
  * `git clone --recurse-submodules https://github.com/bcgsc/btllib` in order to obtain all the code.
  * In `btllib` dir, run `meson build` to create a build directory.
- Every time you want to run tests, in the `build` dir:
  * `ninja wrap` to regenerate wrappers.
  * `ninja test` to build wrappers and tests, and run tests.
- Before making a pull request, in the `build` dir:
  * `ninja quality-assurance` to make sure all CI tests pass.
  * Make a commit after the above step, in case it has made any changes to wrappers or formatting. Don't commit the changes made to the `sdsl-lite` subproject. Meson config file adjusts the `sdsl-lite` config in order for it to work for `btllib`, but this is done ad hoc and is not necessary to be committed. By doing it ad hoc we keep a list of differences compared to the upstream repository.
- Before making a release, in the `build` dir:
  * Do the same as for a pull request and
  * `ninja docs` to regenerate docs to reflect the release and then commit the changes.
  * `meson dist --allow-dirty` to generate a self-contained package based on the last commit. `--allow-dirty` permits making a distributable with uncommited changes. This is necessary as `sdsl-lite` dependency has ad hoc changes made during the build process. The resulting distributable will be compressed with xz. For easier use, decompress it and then compress with gzip. Attach the resulting file to the release.

The following are all the available `ninja` commands which can be run within `build` directory:
- `ninja format` formats the whitespace in code (requires clang-format 8+).
- `ninja wrap` wraps C++ code for Python (requires SWIG 4.0+).
- `ninja clang-tidy` runs clang-tidy on C++ code and makes sure it passes (requires clang-tidy 8+).
- `ninja` builds the tests and wrapper libraries / makes sure they compile.
- `ninja test` runs the tests.
- `ninja sanitize-undefined` runs undefined sanitization.
- `ninja test-wrappers` tests whether wrappers work.
- `ninja docs` generates code documentation from comments (requires Doxygen).
- `ninja quality-assurance` runs `clang-format`, `wrap`, `clang-tidy`, `test`, `sanitize-undefined`, and `test-wrappers`. These are all checked at the CI test.

Credits
---
- Author: [Vladimir Nikolic](https://github.com/vlad0x00)
- Components:
  - [Hamid Mohamadi](https://github.com/mohamadi) and [Parham Kazemi](https://github.com/parham-k) for [ntHash](https://github.com/bcgsc/ntHash)
  - [Justin Chu](https://github.com/JustinChu) for [MIBloomFilter](https://github.com/bcgsc/btl_bloomfilter)
  - [Chase Geigle](https://github.com/skystrife) for [cpptoml](https://github.com/skystrife/cpptoml)
  - Simon Gog, Timo Beller, Alistair Moffat, and Matthias Petri for [sdsl-lite](https://github.com/simongog/sdsl-lite)
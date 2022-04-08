# Attempts to find libc++ and an appropriate ABI (libc++abi or libcxxrt)
# when using clang and libc++ together.

message("-- Locating libc++...")
find_library(LIBCXX_LIBRARY NAMES c++ cxx)
if(LIBCXX_LIBRARY)
  message("-- Located libc++: ${LIBCXX_LIBRARY}")
  set(LIBCXX_OPTIONS "-stdlib=libc++")
  get_filename_component(LIBCXX_LIB_PATH ${LIBCXX_LIBRARY}
    DIRECTORY)
  find_path(LIBCXX_PREFIX c++/v1/algorithm
    PATHS
      ${LIBCXX_LIB_PATH}/../include
      /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include
      /Library/Developer/CommandLineTools/usr/include
      ${CMAKE_SYSTEM_PREFIX_PATH})

  if (LIBCXX_PREFIX)
    set(LIBCXX_INCLUDE_DIR ${LIBCXX_PREFIX}/c++/v1/)
    message("-- Located libc++ include path: ${LIBCXX_INCLUDE_DIR}")
  else()
    message("-- Failed to find libc++ include path!")
  endif()

  message("--     Locating libc++'s abi...")
  find_library(LIBCXXABI_LIBRARY NAMES c++abi)
  find_library(LIBCXXRT_LIBRARY NAMES cxxrt)
  if(LIBCXXABI_LIBRARY)
    message("--     Found libc++abi: ${LIBCXXABI_LIBRARY}")
    set(CXXABI_LIBRARY ${LIBCXXABI_LIBRARY})
  elseif(LIBCXXRT_LIBRARY)
    message("--     Found libcxxrt: ${LIBCXXRT_LIBRARY}")
    set(CXXABI_LIBRARY ${LIBCXXRT_LIBRARY})
  else()
    message("--     No abi library found. "
      "Attempting to continue without one...")
    set(CXXABI_LIBRARY "")
  endif()
else()
  message("-- Could not find libc++!")
endif()

macro(set_libcxx_required_flags)
  if (LIBCXX_OPTIONS)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${LIBCXX_OPTIONS}")
  endif()

  if (CXXABI_LIBRARY)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${LIBCXX_LIB_PATH}")
  endif()
  if (LIBCXX_INCLUDE_DIR)
    set(CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES} ${LIBCXX_INCLUDE_DIR}")
  endif()
endmacro()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBCXX DEFAULT_MSG
  LIBCXX_LIBRARY
  LIBCXX_INCLUDE_DIR
  LIBCXX_LIB_PATH
  LIBCXX_OPTIONS
  CXXABI_LIBRARY)

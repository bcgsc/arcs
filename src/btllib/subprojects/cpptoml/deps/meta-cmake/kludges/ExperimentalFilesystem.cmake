# experimental::filesystem is packaged as a separate static library in
# GCC >= 5.3 as -lstdc++fs
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  message("-- Locating libstdc++ filesystem library")
  find_library(STDCXX_FILESYSTEM
    NAMES stdc++fs
    HINTS ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})

  if (STDCXX_FILESYSTEM)
    message("-- Found libstdc++ filesystem library: ${STDCXX_FILESYSTEM}")
    set(STD_EXPERIMENTAL_FILESYSTEM_LIBRARIES ${STDCXX_FILESYSTEM})
    set(STD_EXPERIMENTAL_FILESYSTEM_FOUND TRUE)
  else()
    message("-- Locating libstdc++ filesystem library - not found")
  endif()
endif()

# experimental::filesystem is packaged in the separate static library in
# libc++ >= 3.9 as -lc++experimental
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND LIBCXX_LIBRARY)
  message("-- Locating libc++experimental library")
  find_library(LIBCXX_EXPERIMENTAL
    NAMES c++experimental
    HINTS ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES} ${LIBCXX_LIB_PATH})

  if (LIBCXX_EXPERIMENTAL)
    message("-- Found libc++experimental library: ${LIBCXX_EXPERIMENTAL}")
    set(STD_EXPERIMENTAL_FILESYSTEM_LIBRARIES ${LIBCXX_EXPERIMENTAL})
    set(STD_EXPERIMENTAL_FILESYSTEM_FOUND TRUE)
  else()
    message("-- Locating libc++experimental library - not found")
  endif()
endif()

# experimental::filesystem is available by default for MSVC 2015 and newer
if (MSVC_VERSION AND NOT (MSVC_VERSION LESS 1900))
  set(STD_EXPERIMENTAL_FILESYSTEM_FOUND TRUE)
endif()

if (STD_EXPERIMENTAL_FILESYSTEM_FOUND)
  message("-- Determining experimental filesystem library capabilities")

  if (STD_EXPERIMENTAL_FILESYSTEM_LIBRARIES)
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} ${STD_EXPERIMENTAL_FILESYSTEM_LIBRARIES}")
  endif()
  check_cxx_source_compiles("
  #include <experimental/filesystem>

  int main()
  {
      std::string path = \"/home/meta/tmp/123456\";
      const auto p1 = std::experimental::filesystem::u8path(path);
      std::experimental::filesystem::remove_all(p1);
      const auto str = p1.u8string();
      return 0;
  }" META_HAS_EXPERIMENTAL_FILESYSTEM)

  if (META_HAS_EXPERIMENTAL_FILESYSTEM)
    target_compile_definitions(compiler-kludges INTERFACE
      -DMETA_HAS_EXPERIMENTAL_FILESYSTEM)

    if (STD_EXPERIMENTAL_FILESYSTEM_LIBRARIES)
      target_link_libraries(compiler-kludges INTERFACE
        ${STD_EXPERIMENTAL_FILESYSTEM_LIBRARIES})
    endif()

    # experimental::filesystem::remove_all doesn't recurse properly as of GCC
    # 5.3.
    set(META_REMOVE_ALL_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/meta-filesystem-test)
    file(MAKE_DIRECTORY ${META_REMOVE_ALL_TEST_DIR})
    file(MAKE_DIRECTORY ${META_REMOVE_ALL_TEST_DIR}/subdir)
    file(WRITE ${META_REMOVE_ALL_TEST_DIR}/subdir/file "Just some testing text")

    check_cxx_source_runs("
    #include <experimental/filesystem>

    int main()
    {
        std::experimental::filesystem::remove_all(\"${META_REMOVE_ALL_TEST_DIR}\");
        return 0;
    }" META_HAS_EXPERIMENTAL_FILESYSTEM_REMOVE_ALL)

    if (META_HAS_EXPERIMENTAL_FILESYSTEM_REMOVE_ALL)
      target_compile_definitions(compiler-kludges INTERFACE
        -DMETA_HAS_EXPERIMENTAL_FILESYSTEM_REMOVE_ALL)
    endif()

    file(REMOVE_RECURSE ${META_REMOVE_ALL_TEST_DIR})
    unset(META_REMOVE_ALL_TEST_DIR)
  endif()
endif()

if (NOT META_HAS_EXPERIMENTAL_FILESYSTEM)
  message("-- No suitable experimental filesystem library found")
endif()

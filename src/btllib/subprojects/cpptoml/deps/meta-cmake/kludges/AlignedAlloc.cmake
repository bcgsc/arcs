set(META_FOUND_ALIGNED_ALLOC_IMPL 0)

check_cxx_source_compiles("
#include <cstdlib>

int main()
{
    ::aligned_alloc(64, 128);
    return 0;
}" META_HAS_ALIGNED_ALLOC)

if (META_HAS_ALIGNED_ALLOC)
  set(META_FOUND_ALIGNED_ALLOC_IMPL 1)
  target_compile_definitions(compiler-kludges INTERFACE
    -DMETA_HAS_ALIGNED_ALLOC)
endif()

if (NOT META_FOUND_ALIGNED_ALLOC_IMPL AND UNIX)
  check_cxx_source_compiles("
  #include <cstdlib>

  int main()
  {
      void* ptr;
      ::posix_memalign(&ptr, 64, 128);
      return 0;
  }" META_HAS_POSIX_MEMALIGN)

  if (META_HAS_POSIX_MEMALIGN)
    set(META_FOUND_ALIGNED_ALLOC_IMPL 1)
    target_compile_definitions(compiler-kludges INTERFACE
      -DMETA_HAS_POSIX_MEMALIGN)
  endif()
endif()

if (NOT META_FOUND_ALIGNED_ALLOC_IMPL AND WIN32)
  check_cxx_source_compiles("
  #include <malloc.h>

  int main()
  {
      ::_aligned_malloc(128, 64);
      return 0;
  }" META_HAS_ALIGNED_MALLOC)

  if (META_HAS_ALIGNED_MALLOC)
    set(META_FOUND_ALIGNED_ALLOC_IMPL 1)
    target_compile_definitions(compiler-kludges INTERFACE
      -DMETA_HAS_ALIGNED_MALLOC)
  endif()
endif()

if (NOT META_FOUND_ALIGNED_ALLOC_IMPL)
  message(FATAL_ERROR "Failed to find a suitable aligned allocation routine")
endif()

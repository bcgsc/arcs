check_cxx_source_compiles("
#include <cstdint>
int main()
{
    uint64_t val = 1;
    auto v = __builtin_ctzll(val);
    return 0;
}" META_HAS_BUILTIN_CTZLL)

if (META_HAS_BUILTIN_CTZLL)
  target_compile_definitions(compiler-kludges INTERFACE
    -DMETA_HAS_BUILTIN_CTZLL)
endif()

check_cxx_source_compiles("
#include <cstdint>
#include <intrin.h>

int main()
{
    uint64_t val = 1;
    unsigned long result;
    if (_BitScanForward64(&result, val))
        return 1;
    return 0;
}" META_HAS_BITSCANFORWARD64)

if (META_HAS_BITSCANFORWARD64)
  target_compile_definitions(compiler-kludges INTERFACE
    -DMETA_HAS_BITSCANFORWARD64)
endif()

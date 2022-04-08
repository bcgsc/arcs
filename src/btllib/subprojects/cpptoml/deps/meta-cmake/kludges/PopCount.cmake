check_cxx_source_compiles("
#include <cstdint>
int main()
{
    uint64_t val = 1;
    auto v = __builtin_popcountll(val);
    return 0;
}" META_HAS_BUILTIN_POPCOUNTLL)

if (META_HAS_BUILTIN_POPCOUNTLL)
  target_compile_definitions(compiler-kludges INTERFACE
    -DMETA_HAS_BUILTIN_POPCOUNTLL)
endif()

check_cxx_source_compiles("
#include <cstdint>
#include <intrin.h>

int main()
{
    uint64_t val = 1;
    auto result = __popcnt64(val);
    return 0;
}" META_HAS_POPCOUNT64)

if (META_HAS_POPCOUNT64)
  target_compile_definitions(compiler-kludges INTERFACE
    -DMETA_HAS_POPCOUNT64)
endif()

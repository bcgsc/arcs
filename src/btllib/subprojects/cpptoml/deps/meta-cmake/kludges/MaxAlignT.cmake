check_cxx_source_compiles("
#include <cstddef>

int main() {
    auto i = alignof(std::max_align_t);
    return 0;
}" META_HAS_STD_MAX_ALIGN_T)

if (META_HAS_STD_MAX_ALIGN_T)
  set(META_MAX_ALIGN_T std::max_align_t)
else()
  set(META_MAX_ALIGN_T ::max_align_t)
endif()

target_compile_definitions(compiler-kludges INTERFACE
                           -DMETA_MAX_ALIGN_T=${META_MAX_ALIGN_T})

check_cxx_source_compiles("
#include <experimental/optional>
int main() {
    std::experimental::optional<int> x;
    return 0;
}" META_HAS_EXPERIMENTAL_OPTIONAL)

if (META_HAS_EXPERIMENTAL_OPTIONAL)
  target_compile_definitions(compiler-kludges INTERFACE
                             -DMETA_HAS_EXPERIMENTAL_OPTIONAL)
endif()

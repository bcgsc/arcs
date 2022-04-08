check_cxx_source_compiles("
#include <memory>
int main() {
    auto i = std::make_unique<int>(1);
    return 0;
}" META_HAS_STD_MAKE_UNIQUE)

if(META_HAS_STD_MAKE_UNIQUE)
  target_compile_definitions(compiler-kludges INTERFACE
                             -DMETA_HAS_STD_MAKE_UNIQUE)
endif()

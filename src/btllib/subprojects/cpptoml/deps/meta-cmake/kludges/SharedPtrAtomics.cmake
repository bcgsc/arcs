check_cxx_source_compiles("
#include <atomic>
#include <memory>
int main () {
    auto sp = std::make_shared<int>(1);
    auto sp2 = std::atomic_load(&sp);
    return 0;
}" META_HAS_STD_SHARED_PTR_ATOMICS)

if(META_HAS_STD_SHARED_PTR_ATOMICS)
  target_compile_definitions(compiler-kludges INTERFACE
                             -DMETA_HAS_STD_SHARED_PTR_ATOMICS=1)
endif()

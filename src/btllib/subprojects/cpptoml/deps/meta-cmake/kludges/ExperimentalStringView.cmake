check_cxx_source_compiles("
#include <experimental/string_view>
int main() {
    const std::experimental::string_view sv = \"hello world\";
    // test that string_view has to_string() const method
    // Xcode 6.4 appears to have shipped a string_view without it
    auto str = sv.to_string();
    return 0;
}" META_HAS_EXPERIMENTAL_STRING_VIEW)

if (META_HAS_EXPERIMENTAL_STRING_VIEW)
  target_compile_definitions(compiler-kludges INTERFACE
                             -DMETA_HAS_EXPERIMENTAL_STRING_VIEW)
endif()

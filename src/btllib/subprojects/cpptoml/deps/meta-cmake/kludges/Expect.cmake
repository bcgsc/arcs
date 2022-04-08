check_cxx_source_compiles("
int main()
{
    long x = 1;
    if (__builtin_expect(x == 1, 0))
        return 1;
    return 0;
}" META_HAS_BUILTIN_EXPECT)

if (META_HAS_BUILTIN_EXPECT)
  target_compile_definitions(compiler-kludges INTERFACE
    -DMETA_HAS_BUILTIN_EXPECT)
endif()

check_cxx_source_compiles("
#include <future>

class no_def_ctor
{
  public:
    no_def_ctor() = delete;
    no_def_ctor(int x) : x_(x) { }
  private:
    int x_;
};

int main()
{
    std::promise<no_def_ctor> p;
    return 0;
}" META_HAS_PROMISE_WITH_NO_DEFAULT_CTOR)

if (META_HAS_PROMISE_WITH_NO_DEFAULT_CTOR)
  target_compile_definitions(compiler-kludges INTERFACE
                             -DMETA_HAS_PROMISE_WITH_NO_DEFAULT_CTOR)
endif()

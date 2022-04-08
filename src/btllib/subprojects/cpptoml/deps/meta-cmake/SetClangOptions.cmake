# Sets a bunch of configuration options when using Clang
macro(SetClangOptions interface_lib)
  # Enable -Wconversion on clang, since it's not *too* noisy there.
  #
  # As of GCC 5.2.0, there are still too many spurious warnings to bother
  # enabling this there.
  target_compile_options(${interface_lib} INTERFACE "-Wconversion")

  if(CMAKE_GENERATOR STREQUAL "Ninja")
    target_compile_options(${interface_lib} INTERFACE "-fcolor-diagnostics")
  endif()

  if (NOT APPLE AND ENABLE_LIBCXX)
    find_package(LIBCXX REQUIRED)
    set_libcxx_required_flags()
  endif()
endmacro()

project('btllib', 'cpp',
        version : '1.3.0',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++11', 'warning_level=3', 'werror=true' ])

# Dependencies and compiler
# ===========================================================

threads_dep = dependency('threads')
openmp_dep = dependency('openmp', required : false)
deps = [ threads_dep, openmp_dep ]

compiler = meson.get_compiler('cpp')
compiler_id = compiler.get_id()
compiler_version = compiler.version()

if compiler_id == 'gcc' and compiler_version.version_compare('<5.0.0')
  error('GCC ' + compiler_version + ' is unsupported.')
endif
if compiler_id == 'clang' and compiler_version.version_compare('<4.0.0')
  error('Clang ' + compiler_version + ' is unsupported.')
endif

fs = import('fs')

# Configuration
# ===========================================================

conf = configuration_data()
conf.set('PROJECT_VERSION', meson.project_version(), description : 'btllib version. Useful to append to the help prompt of recipes.')
configure_file(input : 'include/btllib/btllib_config.hpp.in', output : 'btllib_config.hpp', configuration : conf)

# Compile and link arguments
# ===========================================================

global_args = []
global_link_args = [ '-ldl' ]

if compiler_id == 'clang'
  global_link_args += [ '-lstdc++', '-lm' ]
endif

sdsl_lib = meson.get_compiler('cpp').find_library('sdsl', required : false, dirs : meson.source_root() + '/external/sdsl-lite/installdir/lib')
global_args += [ '-I' + meson.source_root() + '/external/sdsl-lite/installdir/include' ]
global_link_args += [ '-L' + meson.source_root() + '/external/sdsl-lite/installdir/lib', '-lsdsl', '-ldivsufsort', '-ldivsufsort64' ]

add_global_arguments(global_args, language : 'cpp')
add_global_link_arguments(global_link_args, language : 'cpp')

# Custom run commands
# ===========================================================

rootpath = meson.source_root()

run_target('build-sdsl', command: join_paths(rootpath, 'scripts/build-sdsl'))

clang_format = find_program('clang-format', required : false)
if clang_format.found()
  run_target('format', command: join_paths(rootpath, 'scripts/format'))
endif

swig = find_program('swig', required : false)
if swig.found()
  run_target('wrap', command: join_paths(rootpath, 'scripts/wrap'))
endif

clang_tidy = find_program('clang-tidy', required : false)
if clang_tidy.found()
  run_target('tidycheck', command: [ join_paths(rootpath, 'scripts/tidycheck') ] + global_args + global_link_args)
endif

diff = find_program('diff', required : false)
if diff.found() and clang_format.found()
  run_target('formatcheck', command: join_paths(rootpath, 'scripts/formatcheck'))
endif

cppcheck = find_program('cppcheck', required : false)
if cppcheck.found()
  run_target('cppcheck', command: join_paths(rootpath, 'scripts/cppcheck'))
endif

run_target('sanitize-address', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'address'
])
run_target('sanitize-memory', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'memory'
])
run_target('sanitize-thread', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'thread'
])
run_target('sanitize-undefined', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'undefined'
])

doxygen = find_program('doxygen', required : false)
if doxygen.found()
  run_target('docs', command: join_paths(rootpath, 'scripts/docs'))
endif

run_target('complete', command: join_paths(rootpath, 'scripts/complete'))

run_target('test-wrappers', command: join_paths(rootpath, 'scripts/test-wrappers'))

# The library
# ===========================================================

btllib_include = [ include_directories('include'), include_directories('.') ]
btllib_sources = run_command('scripts/get-files', 'src').stdout().strip().split()

btllib_lib = static_library('btllib',
                           btllib_sources,
                           include_directories : btllib_include,
                           dependencies : deps,
                           install : true,
                           install_dir : 'lib')

btllib_dep = declare_dependency(
  link_with : btllib_lib,
  include_directories : btllib_include
)

subdir('wrappers')
subdir('recipes')
if get_option('buildtype') != 'release'
  subdir('tests')
endif

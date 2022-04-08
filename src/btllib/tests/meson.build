if not samtools.found()
  error('tests require samtools to be installed.')
endif

if not gzip.found()
  error('tests require gzip to be installed.')
endif

if not xz.found()
  error('tests require xz to be installed.')
endif

if not bzip2.found()
  error('tests require bzip2 to be installed.')
endif

if not tar.found()
  error('tests require tar to be installed.')
endif

tests_files = run_command('../scripts/get-tests').stdout().strip().split()

foreach file : tests_files
    t = file.split('.cpp')[0].split('.cxx')[0].split('/')[-1]
    test(t, executable(t + '_test', file,
        include_directories : btllib_include,
        dependencies : deps + [ btllib_dep ])
    )
endforeach
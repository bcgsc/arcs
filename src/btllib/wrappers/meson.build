python3 = find_program('python3', required : false)
if python3.found()
    python3_config = find_program('python3-config', required : false)
    if python3_config.found()
        subdir('python')
    else
        message('python3-config not found, skipping Python wrappers.')
    endif
else
    message('python3 not found, skipping wrappers.')
endif

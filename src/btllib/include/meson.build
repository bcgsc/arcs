header_paths = run_command('../scripts/get-include-files').stdout().strip().split()

# Do not install
#foreach header : header_paths
#  relative = run_command('../scripts/get-include-relative', header).stdout().strip()
#  install_headers(files(header), subdir: relative)
#endforeach
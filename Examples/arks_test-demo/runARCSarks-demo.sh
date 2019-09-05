#!/bin/bash
#LEC2018
echo "Ensure the following are in your PATH: arcs, LINKS, arcs-make"
echo "	NOTE: arcs-make Makefile is found in Examples directory of repo"
echo Running ARKS makefile....
../arcs-make arks draft=test_scaffolds reads=test_reads m=50-6000 a=0.9  

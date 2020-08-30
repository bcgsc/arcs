#!/bin/bash

echo Running ARCS makefile
arcs-make arcs-long draft=test_scaffolds reads=test_reads m=8-10000 s=70 a=0.9 l=3 c=3

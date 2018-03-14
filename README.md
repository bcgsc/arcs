![Logo](https://github.com/bcgsc/arcs/blob/master/arcs-logo.png)

# ARCS

Scaffolding genome sequence assemblies using 10X Genomics GemCode/Chromium data.

### Dependencies
* Boost (tested on 1.61)
* GCC (tested on 4.4.7)
* Autotools (if cloning directly from repository) 
* LINKS (tested on 1.8)

### Compilation:
If cloning directly from the repository run:
```
./autogen.sh
```
To compile ARCS run:
```
./configure && make
```
To install ARCS in a specified directory:
```
./configure --prefix=/ARCS/PATH && make install
```
If your boost library headers are not in your PATH you can specify their location:
```
./configure –-with-boost=/boost/path --prefix=/ARCS/PATH && make install
```

### ARCS+LINKS Pipeline 

There are three steps to the pipeline:

1. Run ARCS to generate a Graphviz Dot file (.gv). Nodes in the graph are the sequences to scaffold, and edges show that there is evidence to suggest nodes are linked based on the data obtained from the GemCode/Chromium reads.

2. Run the python script Examples/makeTSVfile.py to generate a file named XXX.tigpair_checkpoint file from the ARCS graph file. The XXX.tigpair_checkpoint file will be provided to LINKS in step 3.

3. Run LINKS with the XXX.tigpair_checkpoint file as input. To do this, the base name (-b) must be set to the same name as XXX.

When using the `-D`/`--dist_est` ARCS option to estimate gap sizes, the user is recommended to use LINKS v1.8.6 or later.

An example bash script on how to run the ARCS+LINKS pipeline can be found at: Examples/pipeline_example.sh

you can test your installation by following instructions at: Examples/arcs_test-demo/README.txt
and compare your output to the files provided at: Examples/arcs_test-demo/output/ 

### Citing ARCS

Paper :
https://doi.org/10.1101/100750


LINKS :
http://www.bcgsc.ca/platform/bioinfo/software/links
https://github.com/warrenlr/LINKS


### License  

ARCS Copyright (c) 2016 British Columbia Cancer Agency Branch.  All rights reserved.

ARCS is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact Patrick Rebstein <prebstein@bccancer.bc.ca>

![Logo](https://github.com/bcgsc/arcs/blob/master/arcs-logo.png)

# ARCS

Scaffolding genome sequence assemblies using 10X Genomics GemCode/Chromium data. 
ARCS can be run in 2 modes:
* [ARCS](https://doi.org/10.1101/100750) (default) uses alignments of linked reads to the input contigs 
* [ARKS](https://doi.org/10.1186/s12859-018-2243-x) (`--arks`) uses exact k-mer mapping to associate linked reads to input contigs

Because ARKS is not dependent on read alignments, it is generally much faster than ARCS. However, ARCS is recommended for use with very fragmented assemblies and/or large genomes.

### Dependencies
* Boost (tested on 1.61)
* GCC (tested on 4.4.7)
* Autotools (if cloning directly from repository) 
* LINKS (tested on 1.8)
* Google SparseHash

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

The ARCS+LINKS pipeline requires two input files:
* Draft assembly fasta file
* Interleaved linked reads file (Barcode sequence expected in the BX tag of the read header or in the form "@readname_barcode" ; Run [Long Ranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) on raw chromium reads to produce this interleaved file)

The Makefile located here: Examples/arcs-make will run the full ARCS pipeline. It will also optionally run the misassembly corrector [Tigmint](https://github.com/bcgsc/tigmint) prior to scaffolding with ARCS.

There are three steps to the pipeline:

1. Run ARCS to generate a Graphviz Dot file (.gv). Nodes in the graph are the sequences to scaffold, and edges show that there is evidence to suggest nodes are linked based on the data obtained from the GemCode/Chromium reads.

2. Run the python script Examples/makeTSVfile.py to generate a file named XXX.tigpair_checkpoint file from the ARCS graph file. The XXX.tigpair_checkpoint file will be provided to LINKS in step 3.

3. Run LINKS with the XXX.tigpair_checkpoint file as input. To do this, the base name (-b) must be set to the same name as XXX.

When using the `-D`/`--dist_est` ARCS option to estimate gap sizes, the user is recommended to use LINKS v1.8.6 or later.

An example bash script on how to run the ARCS+LINKS pipeline can be found at: Examples/pipeline_example.sh

### Running ARCS in default mode

The default mode uses alignments of linked reads to contigs to scaffold the input contigs.

To run the pipeline in default mode, run `Examples/arcs-make arcs`. For more info check `Examples/arcs-make help`.

To run the `arcs` executable in default mode, run `arcs <alignments>`. For descriptions of all arguments, run `arcs --help`.

### Running ARCS in '--arks' mode

To run the pipeline in ARKS mode, run `Examples/arcs-make arks`. For more info check `Examples/arcs-make help`.

To run the `arcs` executable in ARKS mode, run `arcs --arks`. For descriptions of all arguments, run `arcs --help`.

## Demo

You can test your installation by running one of our supplied demos:
* ARCS: `Examples/arcs_test-demo`
* ARKS: `Examples/arks_test-demo`

For both, you can compare your output to the files provided in the `output` folders within the above directories.

### Citing ARKS

Paper :
https://doi.org/10.1186/s12859-018-2243-x

### Citing ARCS

Paper :
https://doi.org/10.1101/100750

**NOTE: The supplementary data and scripts have been moved to http://www.bcgsc.ca/downloads/supplementary/ARCS/**

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

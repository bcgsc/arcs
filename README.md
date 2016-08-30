# ARCS

Scaffolding genome sequence assemblies using 10X Genomics GemCode/Chromium data.

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

### Pipeline 

There are three steps to the pipeline:

1. Run ARCS to generate a Graphviz Dot file (.gv). Nodes in the graph are the sequences to scaffold, and edges in the graph show that there is evidence to suggest nodes are linked based on the data obtained from the GemCode/Chromium reads.

2. Run the python script Examples/makeTSVfile.py to generate a file named XXX.tigpair_checkpoint file from the ARCS graph file. The XXX.tigpair_checkpoint file will be provided to LINKS in step 3.

3. Run LINKS with the XXX.tigpair_checkpoint file as input. To do this, the base name (-b) must be set to the same name as XXX.

An example bash script on how to run the ARCS+LINKS pipeline can be found at: Examples/pipeline_example.sh

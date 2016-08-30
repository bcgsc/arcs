## ARCS+LINKS pipeline example
## ARCS and LINKS must be in your path

if [ $# -ne 2 ]; then
	echo "Usage: $(basename $0) <f> <a>" >&2
    echo "-f  Assembled Sequences to further scaffold (Multi-Fasta format)"
    echo "    NOTE: Sequences must include a unique number (id) in the header"
    echo "-a  File of File Names listing all input BAM alignment files."
    echo "    NOTE: alignments must be sorted in order of name"
    echo "           index must be included in read name e.g read1_indexA"
	exit 1
fi

f=$1; shift
a=$1; shift

## Run ARCS w default params
arcs -f $f -a $a -s 98 -c 5 -l 0 -z 500 -m 50-1000 -d 0 -e 30000 -r 0.05 -i 16 -v 1 

## Run python script makeTSVfile.py to convert ARCS graph output to LINKS XXX.tigpair_checkpoint file format
##  NOTE: XXX must be the same as the base name (-b) for LINKS
graph=$f.scaff_s98_c5_l0_d0_e30000_r0.05_original.gv
python makeTSVfile.py $graph $f.c5_e30000_r0.05.tigpair_checkpoint.tsv $f

## Run LINKS with generated XXX.tigpair_checkpoint file as input
touch empty.fof
LINKS -f $f -s empty.fof -k 20 -b $f.c5_e30000_r0.05 -l 5 -t 2 -a 0.3

#!/bin/bash
#RLW2016
#$ -S /bin/bash
#$ -N ARCS
#$ -q all.q
#$ -pe ncpus 1
#$ -l mem_token=45G,mem_free=45G,h_vmem=45G,excl=true
#$ -j y
#$ -o log/$JOB_NAME.$JOB_ID

if [ $# -ne 4 ]; then
	echo "Usage: $(basename $0) <c> <r> <e> <draft assembly>" >&2
	exit 1
fi

c=$1; shift
r=$1; shift
e=$1; shift
draft=$1;shift

#loose dependency on -i
/usr/bin/time -v -o ARCStime_$c-$r-$e.txt arcs -f $draft -a alignments.fof -s 98 -c $c -l 0 -d 0 -r $r -e $e -v 1 -m 20-10000

#!/bin/bash
#RLW2016

c=$1; shift
draft=$1;shift

l_values=(5)
a_values=(0.9)

for l in "${l_values[@]}"
do
    for a in "${a_values[@]}"
    do
        /usr/bin/time -v -o LINKStime_$c-$l-$a.txt LINKS -f $draft -s empty.fof -k 20 -b links_$c-l$l-a$a -l $l -t 2 -a $a -x 1 #-z 1000 NOTE: THE BASELINE UNDER -b MUST MATCH THAT OF THE FILE CREATED BY makeTSVfile.py 
	#/projects/ABySS/assemblies/EC/SRP000220/abyss-1.5.2/bin/abyss-fac -e 2000YOURGENOMESIZE0000000 -jt 1000 links_$c-l$l-a$a.scaffolds.fa' >> assembly-stats_$c'.md'
    done
done

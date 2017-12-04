#!/usr/bin/env python
## Create a tigpair_checkpoint file from ARCS output that LINKS can utilize
import sys
import re
import argparse

index2scaff_name = {}
links_numbering = {}

def readGraphFile(infile):
    with open(infile, 'r') as f:
        for line in f:
            test = re.match(r"(\d+)\s+\[id=\"?([^\]\"]+)\"?\]", line.rstrip())
            if test:
                index = test.group(1)
                scaff_name = test.group(2)
                if scaff_name not in links_numbering:
                    index2scaff_name[index] = scaff_name

def makeLinksNumbering(scaffolds_fasta):
    counter = 0
    with open(scaffolds_fasta, 'r') as f:
        for line in f:
            if line[0] == ">":
                seq_id = line.rstrip().split()[0][1:]
                counter += 1
                links_numbering[seq_id] = str(counter)


def writeTSVFile(infile, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as w:
            for line in f:
                test = re.search(r"(\d+)--(\d+)\s+\[label=(\d+), weight=(\d+)", line.rstrip())
                if test:
                    scaffA = index2scaff_name[test.group(1)]
                    scaffB = index2scaff_name[test.group(2)]
                    label = int(test.group(3))
                    links = int(test.group(4))

                    if scaffA > scaffB:
                        scaffA, scaffB = scaffB, scaffA

                    outScaffA = ""
                    outScaffB = ""
                    rOutScaffA = ""
                    rOutScaffB = ""
                    outScaffA = outScaffB = rOutScaffA = rOutScaffB = ""

                    # HH : rA->B or rB->A
                    if label == 0:
                        outScaffA = "r" + links_numbering[scaffA]
                        outScaffB = "f" + links_numbering[scaffB]

                        rOutScaffB = "r" + links_numbering[scaffB]
                        rOutScaffA = "f" + links_numbering[scaffA]
                    # HT : rA->rB or B->A
                    elif label == 1:
                        outScaffA = "r" + links_numbering[scaffA]
                        outScaffB = "r" + links_numbering[scaffB]

                        rOutScaffB = "f" + links_numbering[scaffB]
                        rOutScaffA = "f" + links_numbering[scaffA]
                    # TH : A->B or rB->rA
                    elif label == 2:
                        outScaffA = "f" + links_numbering[scaffA]
                        outScaffB = "f" + links_numbering[scaffB]

                        rOutScaffB = "r" + links_numbering[scaffB]
                        rOutScaffA = "r" + links_numbering[scaffA]
                    # TT : A->rB or B->rA
                    elif label == 3:
                        outScaffA = "f" + links_numbering[scaffA]
                        outScaffB = "r" + links_numbering[scaffB]

                        rOutScaffB = "f" + links_numbering[scaffB]
                        rOutScaffA = "r" + links_numbering[scaffA]

                    match = re.search(r"d=(\d+)", line.rstrip())
                    if match:
                        dist = int(match.group(1))
                    else:
                        dist = 10

                    if dist < 0:
                        dist_category = -1
                    elif dist < 500:
                        dist_category = 500
                    elif dist < 1000:
                        dist_category = 1000
                    elif dist < 5000:
                        dist_category = 5000
                    else:
                        dist_category = 10000

                    gap = links * dist

                    string = str(dist_category) + "\t"  + outScaffA + "\t" + outScaffB + "\t" + str(links) + "\t" + str(gap) + "\n"
                    w.write(string)

                    stringR = str(dist_category) + "\t" + rOutScaffB + "\t" + rOutScaffA + "\t" + str(links) + "\t" + str(gap) + "\n"
                    w.write(stringR)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a XXX.tigpair_checkpoint.tsv file from ARCS graph output that LINKS can utilize')
    parser.add_argument('graph_file', help='ARCS graph file output (.gv)', type=str)
    parser.add_argument('out_file', help='Output file name. Must be named XXX.tigpair_checkpoint.tsv, where XXX is same as base name (-b) given to LINKS.', type=str)
    parser.add_argument('fasta_file', help='FASTA file of sequences to scaffold', type=str)
    args = parser.parse_args()

    readGraphFile(args.graph_file)
    makeLinksNumbering(args.fasta_file)
    writeTSVFile(args.graph_file, args.out_file)


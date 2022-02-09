#!/usr/bin/env python
"""
Create a tigpair_checkpoint file from ARCS output that LINKS can utilize
"""
import re
import argparse

index2scaff_name = {}
links_numbering = {}

def readGraphFile(infile):
    """Read ARCS graph file output (.gv)"""
    with open(infile, 'r') as f:
        for line in f:
            test = re.match(r"(\d+)\s+\[id=\"?([^\]\"]+)\"?\]", line.rstrip())
            if test:
                index = test.group(1)
                scaff_name = test.group(2)
                if scaff_name not in links_numbering:
                    index2scaff_name[index] = scaff_name

def makeLinksNumbering(scaffolds_fasta):
    """"""
    counter = 0
    with open(scaffolds_fasta, 'r') as f:
        for line in f:
            if line[0] == ">":
                seq_id = line.rstrip().split()[0][1:]
                counter += 1
                links_numbering[seq_id] = str(counter)


def writeTSVFile(infile, outfile):
    """Create tigpair_checkpoint file"""
    with open(infile, 'r') as fin:
        with open(outfile, 'w') as fout:
            for line in fin:
                test = re.search(r"(\d+)--(\d+)\s+\[label=(\d+), weight=(\d+)", line.rstrip())
                if test:
                    scaffa = index2scaff_name[test.group(1)]
                    scaffb = index2scaff_name[test.group(2)]
                    label = int(test.group(3))
                    links = int(test.group(4))

                    if scaffa > scaffb:
                        scaffa, scaffb = scaffb, scaffa

                    out_scaffa = ""
                    out_scaffb = ""
                    rout_scaffa = ""
                    rout_scaffb = ""
                    out_scaffa = out_scaffb = rout_scaffa = rout_scaffb = ""

                    # HH : rA->B or rB->A
                    if label == 0:
                        out_scaffa = "r" + links_numbering[scaffa]
                        out_scaffb = "f" + links_numbering[scaffb]

                        rout_scaffb = "r" + links_numbering[scaffb]
                        rout_scaffa = "f" + links_numbering[scaffa]
                    # HT : rA->rB or B->A
                    elif label == 1:
                        out_scaffa = "r" + links_numbering[scaffa]
                        out_scaffb = "r" + links_numbering[scaffb]

                        rout_scaffb = "f" + links_numbering[scaffb]
                        rout_scaffa = "f" + links_numbering[scaffa]
                    # TH : A->B or rB->rA
                    elif label == 2:
                        out_scaffa = "f" + links_numbering[scaffa]
                        out_scaffb = "f" + links_numbering[scaffb]

                        rout_scaffb = "r" + links_numbering[scaffb]
                        rout_scaffa = "r" + links_numbering[scaffa]
                    # TT : A->rB or B->rA
                    elif label == 3:
                        out_scaffa = "f" + links_numbering[scaffa]
                        out_scaffb = "r" + links_numbering[scaffb]

                        rout_scaffb = "f" + links_numbering[scaffb]
                        rout_scaffa = "r" + links_numbering[scaffa]

                    match = re.search(r"d=(\d+)", line.rstrip())
                    if match:
                        dist = int(match.group(1))
                        est_dist = True
                    else:
                        dist = 100
                        est_dist = False

                    if dist < 0:
                        dist_category = -1
                    elif dist == 100 and est_dist is False:
                        dist_category = 10
                    elif dist < 500:
                        dist_category = 500
                    elif dist < 1000:
                        dist_category = 1000
                    elif dist < 5000:
                        dist_category = 5000
                    else:
                        dist_category = 10000


                    gap = links * dist

                    string = str(dist_category) + "\t" + out_scaffa + "\t" + \
                             out_scaffb + "\t" + str(links) + "\t" + str(gap) + "\n"
                    fout.write(string)

                    stringr = str(dist_category) + "\t" + rout_scaffb + "\t" + \
                              rout_scaffa + "\t" + str(links) + "\t" + str(gap) + "\n"
                    fout.write(stringr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a XXX.tigpair_checkpoint.tsv file from ARCS graph output '
                                                 'that LINKS can utilize')
    parser.add_argument('graph_file', help='ARCS graph file output (.gv)', type=str)
    parser.add_argument('out_file', help='Output file name. Must be named XXX.tigpair_checkpoint.tsv, where XXX '
                                         'is same as base name (-b) given to LINKS.', type=str)
    parser.add_argument('fasta_file', help='FASTA file of sequences to scaffold', type=str)
    args = parser.parse_args()

    readGraphFile(args.graph_file)
    makeLinksNumbering(args.fasta_file)
    writeTSVFile(args.graph_file, args.out_file)



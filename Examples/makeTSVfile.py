## Create a tigpair_checkpoint file from ARCS output that LINKS can utilize
import sys
import re
import argparse

index2scaff_name = {} 
links_numbering = {}

def readGraphFile(infile):
    with open(infile, 'r') as f:
        for line in f:
            test = re.match(r"(\d+)\s+\[id=(\d+)\]", line.rstrip())
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
                test = re.match(r">(\d+)", line.rstrip())
                counter += 1
                links_numbering[test.group(1)] = str(counter)
                
                    
def writeTSVFile(infile, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as w:
            for line in f:
                test = re.match(r"(\d+)--(\d+)\s+\[label=(\d+), weight=(\d+)\]", line.rstrip())
                if test:
                    scaffA = index2scaff_name[test.group(1)]
                    scaffB = index2scaff_name[test.group(2)]
                    label = int(test.group(3))
                    links = int(test.group(4))
                    
                    if int(scaffA) > int(scaffB):
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

                    gap = links*10

                    string = str(500) + "\t"  + outScaffA + "\t" + outScaffB + "\t" + str(links) + "\t" + str(gap) + "\n"
                    w.write(string)

                    stringR = str(500) + "\t" + rOutScaffB + "\t" + rOutScaffA + "\t" + str(links) + "\t" + str(gap) + "\n"
                    w.write(stringR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a XXX.tigpair_checkpoint file from ARCS graph output that LINKS can utilize')
    parser.add_argument('graph_file', help='ARCS graph file output (.gv)', type=str)
    parser.add_argument('out_file', help='Output file name. Must be named XXX.tigpair_checkpoint.tsv, where XXX is same as base name (-b) given to LINKS.', type=str)
    parser.add_argument('fasta_file', help='FASTA file of sequences to scaffold', type=str)
    args = parser.parse_args()

    readGraphFile(args.graph_file)
    makeLinksNumbering(args.fasta_file)
    writeTSVFile(args.graph_file, args.out_file)


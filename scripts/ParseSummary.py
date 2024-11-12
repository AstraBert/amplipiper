import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import os

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--path", dest="IN", help="Input path")
parser.add_option("--primer", dest="PR", help="primers file")
parser.add_option("--samples", dest="SA", help="samples file")
parser.add_option("--output", dest="OUT", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def getIDS(x):
    ''' get IDS '''
    L = []
    for l in load_data(x):
        a = l.rstrip().split(",")
        if l.startswith("ID,"):
            continue
        L.append(a[0])
    return L


TEST = {2:
        {2: [(0.5, 0.5)],
         3: [(0.3333333333333333, 0.6666666666666666)],
         4: [(0.5, 0.5), (0.25, 0.75)]},
        3:
        {3: [
            (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
            4: [(0.25, 0.25, 0.5)]},
        4: {4: [(0.25, 0.25, 0.25, 0.25)]}}

# Perform chi-square goodness-of-fit test
# chi2_stat, p_val = chisquare(observed, f_exp=expected)

print("ID,Locus,HaplotypesCount,TotalReads,ReadCount,FrequencyOfTotal,HaplotypesLengths")
for LOCUS in getIDS(options.PR):
    for IND in getIDS(options.SA):
        if not os.path.exists(f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/results.txt"):
            RD, PC, UR, SEQL = [], "NA", "NA", "NA"
            print(",".join([IND,
                            LOCUS,
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "NA"]))
        else:
            FILE = load_data(
                f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/results.txt")
            RD, PC = [], []
            for l in FILE:
                if l.startswith("-->"):
                    RD.append(l.split("contains ")[1].split(" sequences")[0])
                    PC.append(
                        str(round(float(l.split("sequences (")[1].split("% of total")[0])/100, ndigits=3)))
                if l.startswith("- used_reads"):
                    URline = l.rstrip().split("= ")
                    if len(URline) > 1:
                        UR = URline[1]
                    else:
                        UR = "0"
            if os.path.exists(f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/{LOCUS}_consensussequences.fasta"):
                FILE = load_data(
                    f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/{LOCUS}_consensussequences.fasta")
                SEQL = []
                for l in FILE:
                    if l.startswith(">"):
                        continue
                    SEQL.append(str(len(l.rstrip())))
            else:
                SEQL = ["NA"]
            HAP = len(RD)
            if RD == []:
                RD = ["NA"]
            if PC == []:
                PC = ["NA"]

            print(",".join([IND,
                            LOCUS,
                            str(HAP),
                            UR,
                            "/".join(RD),
                            "/".join(PC),
                            "/".join(SEQL)]))

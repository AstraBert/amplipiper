from optparse import OptionParser, OptionGroup
import numpy as np
from scipy.stats import multinomial
from scipy.stats import chisquare
from collections import defaultdict as d
import os

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--path", dest="PA", help="Input file")
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


def likelihood(observed, expected):
    # Multinomial likelihood: multinomial.pmf(observed counts | expected proportions)
    # print(observed, np.sum(observed), expected)
    return multinomial.pmf(observed, n=np.sum(observed), p=expected)


def estimatePloidy(OBS, output, Locus, ID):
    PLO = d(lambda: d(list))
    OBS = [int(x) for x in OBS]
    COUNTl = len(OBS)
    for k, v in TEST[COUNTl].items():
        for EXP in v:
            EXP /= np.sum(EXP)
            LH = likelihood(sorted(OBS), EXP)
            PLO[LH]["samples"] = sorted(OBS)
            PLO[LH]["PLOIDY"].append([str(COUNTl), str(k)])
            output.write(",".join(
                [ID, Locus, "/".join([str(x) for x in OBS]), str(COUNTl), str(k), str(LH)])+"\n")
    # OBSold = OBS
    # while (len(OBS) >= 3):
    #     Total = sum(OBS)
    #     OBSold = sorted(OBSold)[1:]
    #     OBS = sorted(OBS)[1:]
    #     Freq = [x/sum(OBS) for x in OBS]
    #     OBS = [round(x*Total, 0) for x in Freq]
    #     COUNTl = len(OBS)

    #     for k, v in TEST[COUNTl].items():
    #         for EXP in v:
    #             EXP /= np.sum(EXP)
    #             LH = likelihood(sorted(OBS), EXP)
    #             PLO[LH]["samples"].extend(OBSold)
    #             PLO[LH]["PLOIDY"].append([str(COUNTl), str(k)])
    #             output.write(",".join(
    #                 [ID, Locus, "/".join([str(x) for x in OBS]), str(COUNTl), str(k), str(LH)])+"\n")
    return PLO


TEST = {2:
        {1: [np.array([0.0, 1.0])],
         2: [np.array([1/2, 1/2])],
         3: [np.array([1/3, 2/3])],
         4: [np.array([1/2, 1/2]), np.array([1/4, 3/4])]},
        3:
        {1: [np.array([0.0, 0.0, 1.0])],
         2: [np.array([0.0, 1/2, 1/2])],
            3: [np.array([1/3, 1/3, 1/3])],
            4: [np.array([1/4, 1/4, 1/2])]},
        4: {1: [np.array([0.0, 0.0, 0.0, 1.0])],
            2: [np.array([0.0, 0.0, 1/2, 1/2])],
            3: [np.array([0, 1/3, 1/3, 1/3])],
            4: [np.array([1/4, 1/4, 1/4, 1/4])]}}


# Perform chi-square goodness-of-fit test
# chi2_stat, p_val = chisquare(observed, f_exp=expected)

# dictionary to store the consensus IDs to keep for each locus and sample
KEEP = d(lambda: d(list))
PLOIDY = d(lambda: d(list))
output = open(options.IN.split(".csv")[0]+".likelihoods", "wt")
output.write("ID,Locus,Reads,Alleles,Ploidy,Likelihood\n")
for l in load_data(options.IN):
    # skip header
    if l.startswith("ID,"):
        continue
    ID, Locus, HaplotypesCount, TotalReads, ReadCount, FrequencyOfTotal, HaplotypesLengths = l.rstrip().split(",")
    # skip if NA
    if ReadCount == "NA":
        continue
    # keep single haplotype if homozygous
    Reads = ReadCount.split("/")
    READID = dict(zip(*[Reads, range(len(Reads))]))
    if len(Reads) > 4:
        Reads = sorted(Reads)[:4]
    if len(Reads) == 1:
        KEEP[Locus][ID].append(0)
        PLOIDY[Locus][ID].append([1, 1])
        continue
    # calculate total reads across all consensus sequences
    GrandTotal = 0
    for C in Reads:
        GrandTotal += int(C)
    # calculate frequencies and remove consensus sequences with rel. Freq <10%
    Keep = []
    for C in range(len(Reads)):
        RF = int(Reads[C])/GrandTotal
        if RF < 0.1:
            continue
        Keep.append(C)
    if len(Keep) == 1:
        KEEP[Locus][ID].append(Keep[0])
        PLOIDY[Locus][ID].append([1, 1])
        continue

    NEW = estimatePloidy([Reads[x] for x in Keep], output, Locus, ID)
    BEST = max(NEW.keys())

    for sample in list(set(NEW[BEST]["samples"])):
        KEEP[Locus][ID].append(READID[str(sample)])
    for ploidy in NEW[BEST]["PLOIDY"]:
        PLOIDY[Locus][ID].append(ploidy)
# now print the ploidies for each locus+sample
Out = open(options.IN+".ploidy", "wt")

for l in load_data(options.IN):
    # skip header
    if l.startswith("ID,"):
        Out.write(l.rstrip()+",NoAlleles,ExpectedPloidy\n")
        continue
    ID, Locus, HaplotypesCount, TotalReads, ReadCount, FrequencyOfTotal, HaplotypesLengths = l.rstrip().split(",")
    # print(ReadCount)
    if ReadCount == "NA":
        Out.write(l.rstrip()+",NA,NA\n")
        continue
    if ID in PLOIDY[Locus]:
        NoA, EP = list(zip(*PLOIDY[Locus][ID]))
        Out.write(l.rstrip()+","+"/".join(list(set([str(x) for x in NoA]))) +
                  ","+"/".join([str(x) for x in EP])+"\n")
    else:
        Out.write(l.rstrip()+",NA,NA\n")

Out.close()

##
for Locus, v in KEEP.items():
    directory = options.OUT+"/"+Locus+"/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    Out = open(directory+Locus+".fasta", "wt")
    for ID, v1 in v.items():
        FILE = options.PA+"/"+ID+"/"+Locus+"/"+Locus+"_consensussequences.fasta"
        FASTA = d(str)
        C = ""
        for l in load_data(FILE):
            if l.startswith(">"):
                if C == "":
                    C = 0
                else:
                    C += 1
                continue
            FASTA[C] += l
        C = 1
        for key in v1:
            Out.write(">"+ID+"_"+str(C)+"\n")
            Out.write(FASTA[key])
            C += 1
    Out.close()

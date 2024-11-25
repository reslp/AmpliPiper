from module import load_data
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--primername", dest="PI", help="Input file")
parser.add_option("--names", dest="names", help="Output file")
parser.add_option("--outgroup", dest="OG", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

OG = options.OG.split(",")
OGnew = []

TEST = 0
NAME = d(list)
for l in load_data(options.names):

    # skip header
    if l.startswith("SAMPLE"):
        continue

    # split by comma
    a = l.rstrip().rstrip(",").split(",")

    # First test if Filename contains "Freq_" indicating --freqthreshold 0 in this case do not split by different Consensus sequences and do not append Species names to loci other than the barcode loci
    if "Freq_" in a[0]:
        TEST = 1
        ID = a[0]
    else:
        ID = "-".join(a[0].split("-")[:-1])
        # ID = a[0]

    # If no Species was identified, use "NA"
    if len(a) == 1:
        NAME[ID].append("NA")
        continue

    # If at least one Species identified, parse Species table and only keep best Hit
    HASH = d(lambda: d(lambda: d(str)))
    for sample in a[1:]:

        # replace spaces with "_"
        Spec = sample.split(" (")[0].replace(" ", "_")

        # only keep three letters from genus name and connect with dot to species names
        Spec = Spec.split("_")[0][:3]+"."+Spec.split("_")[1]

        # keep similarity values
        Sim = float(sample.split(" (")[1].split("%)")[0])

        # keep count values
        count = int(sample.split("count = ")[1])

        # fill hash and recode similarity in percent
        HASH[ID][Sim][count] = [Spec, str(round(Sim/100, 2))]

    # obtain highest counts and similarites
    Sim = max(HASH[ID].keys())
    count = max(HASH[ID][Sim].keys())

    # for each ID, retain Spec with highest counts/similarity
    NAME[ID].append(HASH[ID][Sim][count])

    # rename outgroup samples
    for outgroup in OG:
        if outgroup not in a[0]:
            continue
        ID = "-".join(a[0].split("-")[:-1])
        EXT = a[0].split("-")[-1]
        SPEC = HASH[ID][Sim][count]
        OGnew.append(ID+"-"+"_".join(SPEC)+"-"+EXT)


# if --freqthreshold==0, only append names to COX1, etc.
if TEST == 1:
    tree = load_data(options.IN).readline()
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        for k, v in NAME.items():
            # split into ID and number of consensus
            ID = "-".join(k.split("-")[:-1])
            EXT = k.split("-")[-1]
            # replace the ID in tree with the extended ID
            tree = tree.replace(ID + "-"+EXT, ID+"-"+"_".join(v[0])+"-"+EXT)

else:
    # make more specific dictionary for each consensus per locus separately
    NAME2 = d(lambda: d(str))
    # if barcoding locus, keep species ID specific for each consensus number
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        for k, v in NAME.items():
            for i in range(len(v)):
                NAME2[k][str(i+1)] = "_".join(v[i])
    else:
        # if there are more than two consensus sequences for the barcoding locus, skip adding species name as a whole
        for k, v in NAME.items():
            Spec, Sim = list(zip(*v))
            if len(list(set(Spec))) > 1:
                NAME2[k]["N"] = "NA"
            else:
                NAME2[k]["N"] = "_".join(v[0])
    # OK, now read the tree file
    tree = load_data(options.IN).readline()

    # If barcoding locus, attach species names for each consensus
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        for k, v in NAME2.items():
            for I, v1 in v.items():
                tree = tree.replace(k+"-"+I, k+"-"+v1+"-"+I)
    else:
        # for all other, attach species name if only single consensus of diagnostic barcoding locus with species names
        for k, v in NAME2.items():
            tree = tree.replace(k, k+"-"+v["N"])

# print new tree
out = open(options.IN, "w")
out.write(tree+"\n")

# print new IDs in outgrouplist
print(",".join(OGnew))

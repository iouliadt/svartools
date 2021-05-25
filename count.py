from collections import defaultdict
from Bio import SeqIO
from math import *
import svartools as svt
from dictionaries import codontable, polarity


# * Open alignment file
input_file = open(svt.args.alignment)

# write to dictionary
fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

# If the user doesn't specify a reference id
# then the first entry of the alignment is used as ref

if svt.args.reference is None:

    print("Warning: you have not specified a reference identifier"
          " to tabulate the mutations against, the coordinates"
          " will be relative to the alignment position.")

    refID = list(fasta_dict.keys())[0]
else:
    refID = svt.args.reference

# Get reference sequence using the refID key
ref = fasta_dict[refID].seq


def log_base(base, value):

    return(log(value)/log(base))


def increment(dict, key):

    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1


# Create reference codons
for i in range(0, len(ref), 3):
    rcodon = str(ref[i]+ref[i+1]+ref[i+2])


# Implementation of perl's autovivification feature
class AutoVivification(dict):

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


# Define dictionaries using AutoVivification()
dnds = AutoVivification()
mutfreq = AutoVivification()
entropy = AutoVivification()
nucfreq = AutoVivification()
NSnucfreq = AutoVivification()
NSsub_mat = AutoVivification()
seqcnt = 0

# Loop through query sequences
for seq in fasta_dict.keys():

    refgap = 0
    id = seq
    query = list(fasta_dict[seq])

# Loop through reference sequence
    for i in range(0, len(ref), 3):

        rcodon = ref[i]+ref[i+1]+ref[i+2]
        qcodon = query[i]+query[i+1]+query[i+2]

        if rcodon not in codontable.keys():
            # Translate codons not in codontable (like Ns) to 'X'
            raa = 'X'
        else:
            raa = codontable[rcodon]
            if (raa == '-'):
                refgap = refgap + 3

        if qcodon not in codontable.keys():
            #  print("No translation for codon {}".format(qcodon))
            qaa = 'X'
        else:
            qaa = codontable[qcodon]

        aapos = int(((i-refgap)/3)+1)
        mut = str(raa)+str(aapos)+str(qaa)

        if (raa != '-'):
            increment(mutfreq[aapos], mut)

        if raa == qaa and rcodon != qcodon and raa != '-':
            # if the codons are not equal then store the changes e.g. A=>T as synonymous

            for cpos in range(3):

                if ref[cpos+i] != query[cpos+i]:
                    sub = str(ref[cpos+i])+"=>"+str(query[cpos+i])
                    increment(NSnucfreq[sub], "syn")
                    increment(nucfreq[cpos+i][aapos][query[cpos+i]], "syn")

                if ref[cpos+i] == query[cpos+i]:
                    increment(nucfreq[cpos+i][aapos][query[cpos+i]], "syn")

            increment(dnds[aapos], "syn")

            for j in range(3):

                if ref[i+j] != query[i+j] and raa != '-' or '*' and qaa != '-' or '*':
                    sub = str(ref[i+j]+query[i+j])

                    increment(NSsub_mat[sub]["syn"], aapos)

        # if the codons are not equal then store the changes e.g. A=>T as NON-synonymous
        elif raa != qaa and raa != '-':

            for cpos in range(3):

                if ref[cpos+i] != query[cpos+i]:
                    sub = ref[cpos+i]+"=>"+query[cpos+i]

                    increment(NSnucfreq[sub], "nonsyn")
                    increment(nucfreq[cpos+i][aapos][query[cpos+i]], "nonsyn")

                if ref[cpos+i] == query[cpos+i]:
                    increment(nucfreq[cpos+i][aapos][query[cpos+i]], "nonsyn")

            increment(dnds[aapos], "nonsyn")

            for j in range(3):

                if ref[i+j] != query[i+j] and raa != '-' or '*' and qaa != '-' or '*':
                    sub = str(ref[i+j])+str(query[i+j])
                    increment(NSsub_mat[sub]["nonsyn"], aapos)

        elif rcodon == qcodon:
            # if codons are equal we just need to count the freq of each nuc
            for cpos in range(3):
                if ref[cpos+i] == query[cpos+i]:
                    increment(nucfreq[cpos+i][aapos][query[cpos+i]], "syn")

    seqcnt += 1


# Create file name according to user input
if svt.args.outfile:
    fileName = svt.args.outfile
else:
    fileName = svt.args.alignment

# Open files and write headers
if svt.args.aa:
    with open("{}_entropy.txt".format(fileName), 'a') as ENT, open("{}_mutfreq.txt".format(fileName), 'a') as FREQ:

        ENT.write("aa_position\tRef_aa\tentropy\tentropy(base 2)\tnb_nonsyn\tnb_syn\tdN/dS\n")

        FREQ.write("Mutation\tRef_aa\tRef_polarity\taa_position\tQuery_aa\tQuery_polarity\tfrequency\tcount\n")

elif svt.args.nuc:
    with open("{}_nucfreq.txt".format(fileName), 'a') as NUC:

        NUC.write("NucSite\taa_position\tRefNuc\tAsyn\tAnonsyn\tTsyn\tTnonsyn \tGsyn\tGnonsyn\tCsyn\tCnonsyn\tentropy\tRef_polarity\n")


ent = AutoVivification()
ref_polarity = AutoVivification()
shannon = ''
shannon2 = ''

for aapos in sorted(mutfreq.keys()):

    for mut in mutfreq[aapos].keys():

        p = mutfreq[aapos][mut]/seqcnt
        shannon = -p*log(p)/seqcnt

        if shannon == -0.0:
            shannon = 0

        ent[aapos] = shannon
        shannon2 = -(log_base(2, p))*p

        if shannon2 == -0.0:
            shannon2 = 0

        residue = list(mut.split(str(aapos)))
        raa = residue[0]
        ref_polarity[aapos] = polarity[residue[0]]

        if residue[0] != residue[1]:
            if svt.args.aa:
                with open("{}_mutfreq.txt".format(fileName), 'a') as FREQ:
                    FREQ.write(mut + "\t" + residue[0] + "\t" + polarity[residue[0]]
                               + "\t" + str(aapos) + "\t" + residue[1] + "\t" + polarity[residue[1]]
                               + "\t" + str(p) + "\t" + str(mutfreq[aapos][mut]) + "\n")

                if isinstance(dnds[aapos]["syn"], AutoVivification):
                    dnds[aapos]["syn"] = 0
                if isinstance(dnds[aapos]["nonsyn"], AutoVivification):
                    dnds[aapos]["nonsyn"] = 0
                with open("{}_entropy.txt".format(fileName), 'a') as ENT:
                    if dnds[aapos]["syn"] > 0:

                        ENT.write(str(aapos) + '\t' + raa + '\t' + str(shannon) + '\t' + str(shannon2) + '\t' +
                                  str(dnds[aapos]["nonsyn"]) + "\t" + str(dnds[aapos]["syn"]) +
                                  "\t" + str(dnds[aapos]["nonsyn"]/dnds[aapos]["syn"]) + "\n")

                    else:
                        ENT.write(str(aapos) + '\t' + raa + '\t' + str(shannon) + '\t' + str(shannon2) + '\t' +
                                  str(dnds[aapos]["nonsyn"]) + "\t" + str(dnds[aapos]["syn"]) + "\t" +
                                  "inf" + "\n")

if svt.args.nuc:

    sub_mat = AutoVivification()
    totalnucs = 0

    with open("{}_nucfreq.txt".format(fileName), 'a') as NUC:

        NUC.write("NucSite\taa_position\tRefNuc\tAsyn\tAnonsyn\tTsyn\tTnonsyn\tGsyn\tGnonsyn\tCsyn\tCnonsyn\tentropy\tRef_polarity\n")

        for i in range(len(ref)):
            for aapos in nucfreq[i].keys():
                site = i + 1
                if nucfreq[i][aapos]["A"]["nonsyn"] == {}:
                    nucfreq[i][aapos]["A"]["nonsyn"] = 0
                if nucfreq[i][aapos]["T"]["nonsyn"] == {}:
                    nucfreq[i][aapos]["T"]["nonsyn"] = 0
                if nucfreq[i][aapos]["G"]["nonsyn"] == {}:
                    nucfreq[i][aapos]["G"]["nonsyn"] = 0
                if nucfreq[i][aapos]["C"]["nonsyn"] == {}:
                    nucfreq[i][aapos]["C"]["nonsyn"] = 0
                if nucfreq[i][aapos]["A"]["syn"] == {}:
                    nucfreq[i][aapos]["A"]["syn"] = 0
                if nucfreq[i][aapos]["T"]["syn"] == {}:
                    nucfreq[i][aapos]["T"]["syn"] = 0
                if nucfreq[i][aapos]["G"]["syn"] == {}:
                    nucfreq[i][aapos]["G"]["syn"] = 0
                if nucfreq[i][aapos]["C"]["syn"] == {}:
                    nucfreq[i][aapos]["C"]["syn"] = 0

            NUC.write(str(site) + "\t" + str(aapos) + "\t" + ref[i] + "\t" +
                      str(nucfreq[i][aapos]["A"]["syn"]) + "\t" +
                      str(nucfreq[i][aapos]["A"]["nonsyn"]) + "\t")

            NUC.write(str(nucfreq[i][aapos]["T"]["syn"]) + "\t" +
                      str(nucfreq[i][aapos]["T"]["nonsyn"]) + "\t")

            NUC.write(str(nucfreq[i][aapos]["G"]["syn"]) + "\t" +
                      str(nucfreq[i][aapos]["G"]["nonsyn"]) + "\t")

            NUC.write(str(nucfreq[i][aapos]["C"]["syn"]) + "\t" +
                      str(nucfreq[i][aapos]["C"]["nonsyn"]) + "\t")

            NUC.write(str(ent[aapos]) + "\t" + str(ref_polarity[aapos]) + "\n")

            for nuc in nucfreq[i][aapos].keys():
                if nucfreq[i][aapos][nuc] != {}:
                    increment(sub_mat[ref[i]], nuc)

        nucs = nucfreq[i].keys()


# Total number of sites where particular substitutions are observed
    with open("{}_nucfreq.sum.txt".format(fileName), "a") as SUMNUC:
        SUMNUC.write("RefNuc\tdbNuc\tNbSitesWithChange\n")

        for refnuc in sub_mat.keys():

            for dbnuc in sub_mat[refnuc].keys():

                sub = str(refnuc + dbnuc)
                nbnonsyn = NSsub_mat[sub]["nonsyn"].keys()
                nbnonsyn = NSsub_mat[sub]["syn"].keys()
                SUMNUC.write(str(refnuc) + "\t" + str(dbnuc) + "\t" +
                             str(sub_mat[refnuc][dbnuc]) + "\n")

if svt.args.nuc:

    # Total number of mutations that are in non-synonymous or synonymous codons
    with open("{}_dnds.sum.txt".format(fileName), "a") as SUMDNDS:
        SUMDNDS.write("RefNuc\tdbNuc\tNonSyn\tSyn\n")

        for sub in NSnucfreq.keys():

            if "-" not in sub:

                refnuc = sub.split("=>")[0]
                dbnuc = sub.split("=>")[1]

                if NSnucfreq[sub]["nonsyn"]:
                    ns_cnt = NSnucfreq[sub]["nonsyn"]
                else:
                    ns_cnt = 0
                if NSnucfreq[sub]["syn"]:
                    s_cnt = NSnucfreq[sub]["syn"]
                else:
                    s_cnt = 0

                SUMDNDS.write(refnuc + "\t" + dbnuc +
                            "\t" + str(ns_cnt) + "\t" + str(s_cnt) + "\n")

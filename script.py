#!/usr/bin/env python3
"""
Description of the script

The fasta file used for this script was directly download from Arabidopsis.org in April 2023
The input was CDS of arabidopsis that can be found on the website download/sequences section,
URL: https://www.arabidopsis.org/download_files/Sequences/Araport11_blastsets/Araport11_cds_20220914.gz
This script has not been tested with other types of inputs.
"""

import argparse
import re

parser = argparse.ArgumentParser(
    prog="help",
    description="PCR primer generator script that gets a fasta file as an input, \n"
                "and returns forward and reverse primers along with their Tm and GC content \n"
                "You have to provide a fasta file containing CDS sequences downloaded from arabidopsis website \n"
                "E.g. URL: https://www.arabidopsis.org/download_files/Sequences/Araport11_blastsets/Araport11_cds_20220914.gz \n"
                "Don't forget to unzip the file\n",
    epilog="If there are problems contact amirreza.alise@gmail.com"
)
parser.add_argument("-f", "--FastaFile", required=True,
                    help="Input the path for fasta file (MUST contains CDS)\n"
                         "The file must be a text format, either .txt or .fa/.fasta or similar \n"
                         "E.g <folder path>/Araport11_cds_20220914.txt")
parser.add_argument("-k", "--KeyWord", required=True,
                    help="The keyword(s) that you want to search for in the genes. \n"
                         "remember to separate your keywords with comma (,) and no spaces if there are more than one\n"
                         "E.g. SCARECROW,scarecrow")
args = parser.parse_args()

# Opening the file which is given by the user
file = open(args.FastaFile, "r")
sequences = file.read()

file.close()

# if the file doesn't start with a > it means it's not a fasta file!
if not sequences.startswith('>'):
    raise TypeError("Not a Fasta, Check file and start again.")

# Defining the regex to find each line including details of the sequence with given keyword(s)
try: # tries to see if there are more than one keyword
    keyword_list = str(args.KeyWord).split(",") # splits the keywords into a list
    keyword = "|".join(keyword_list) # makes a new string and separates keywords with |
except: # if there is only one keyword as an input
    keyword = str(args.KeyWord)

regex = r'>.*(%s).*\d$' % keyword # %s will be replaced with keyword in regex
list1 = re.finditer(regex, sequences, re.MULTILINE)  ## list of all lines with gene details that has the word SCARECROW

count = 0
pos_start_list = []
pos_end_list = []
genes_list = []

print("These sequences where found with Keyword: ", keyword)
# The following commands extract the pure(cleaned) sequence(Only bases) from the file
for ind_s, i in enumerate(list1):
    string = i.group()
    genes_list.append(string)
    print(ind_s, '.', string)
    try:
        start = int(i.start())
        end = int(i.end())
        pos_start_list.append(start)
        pos_end_list.append(end)
        count += 1
    except:
        print("This is here to keep the flow of the script, just ignore it")

cDNA_list = []
regex_seq = r'(?<=\n)\w.[^>]*(?=>)'  # Starts with a word(ACTG) preceded by a new line and followed by '>'

for j in range(count):
    start1 = pos_end_list[j]  # start of the gene
    seq_temp = sequences[start1 - 2:]

    # starts searching from  the end of keyword(SCARECROW) sequence
    seq_dirty = re.findall(regex_seq, seq_temp)[0]  # returns the first match
    seq_clean = re.sub(r'\n', '', seq_dirty)
    cDNA_list.append(seq_clean)

# I have extracted the sequences and saved them in a list argument (cDNA_list)

'''
Aspects of a primer:
.Length 20-23

.Annealing and Melting Temperatures
Annealing temperature should be 5 degrees below  melting temperature
Forward and Backward strand should have similar melting temperatures

.GC content

.secondary structure

the anealing region must be unique

'''


# Functions to calculate Tm and GC
def tm_gc(primer_seq):
    Tm = 0
    GC_content = 0

    try:
        # Counting each base in primer sequence
        wA = primer_seq.count("A")
        xT = primer_seq.count("T")
        yG = primer_seq.count("G")
        zC = primer_seq.count("C")

        values = (wA, xT, yG, zC)
        # ignoring 0 values because some sequences might have ATG at the begining and it causes the program to return 0
        if all([v != 0 for v in values]):

            # Calculating Melting Tempreture (Tm):
            if len(primer_seq) < 14:
                # TODO For sequences less than 14 nucleotides the formula is
                Tm = (wA + xT) * 2 + (yG + zC) * 4

            if len(primer_seq) >= 14:
                # TODO For sequences longer than 14 nucleotides, the equation used is
                Tm = 64.9 + 41 * (yG + zC - 16.4) / (wA + xT + yG + zC)

            # calculating the GC content of the sequence:
            GC_content = ((yG + zC) / len(seq)) * 100
    except:
        raise "Unexpected Error!"

    return Tm, GC_content


# function to make complement sequence
def complement(sequence):
    # making complement of the string (Forward primer sequence)
    temp = sequence.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    complement_seq = temp.upper()

    return complement_seq


def replacenth(str1, sub, wanted, n):
    """
    :param str1: 'ababababababababab'
    :param sub: 'ab'
    :param wanted: 'CD'
    :param n: 5
    :return: ababababCDabababab
    """
    where = [m.start() for m in re.finditer(sub, str1)][n - 1]
    before = str1[:where]
    after = str1[where:]
    after = after.replace(sub, wanted, 1)
    newString = before + after
    return newString


# we try to find the best primer candidates
# we try with each length (20,21,22,23)
lengths = [20, 21, 22, 23]
# accessing each sequence
counter = 0
for index, seq in enumerate(cDNA_list):
    # trying each length
    f_primer_list = []
    r_primer_list = []
    for l in lengths:
        # searching for sequence with length l after ATG
        # Finding the start codon
        regex_start = r'ATG(.){%s}' % l
        f_primer_seq = re.search(regex_start, seq).group()
        # Making complement sequence (Primer)
        f_primer = complement(f_primer_seq)
        # Adding it to a list of candidates
        f_primer_list.append(f_primer)

        # f_primer_pattern = re.compile(regex_start)
        # ind_ATG = f_primer_pattern.search(seq).start()

        regex_stop = r'.{%s}(?<=(TAG|TGA|TAA))$' % l
        r_primer_seq = re.search(regex_stop, seq).group()

        # Reverse Sequence
        r_primer_rev = r_primer_seq[::-1]

        # Complement of Reversed sequence
        r_primer = complement(r_primer_rev)
        r_primer_list.append(r_primer)

    cand = 0
    for r in r_primer_list:

        # calculating GC content and Tm of the primer sequence
        r_Tm, r_GC = tm_gc(r)

        for f in f_primer_list:
            # calculating GC content and Tm of the primer sequence
            f_Tm, f_GC = tm_gc(f)

            # Checking if the Tm is in the desired range (55 to 62) and not diverging more than 4 degrees from each other
            if 55 <= f_Tm <= 62 and 55 <= r_Tm <= 62 and abs(f_Tm - r_Tm) <= 4:
                cand += 1
                counter += 1
                print(
                    f"____________________ {counter}. PRIMERS FOUND MATCHING THE GIVEN CRITERIA _____________________\n"
                    f"{cand}th candidate primer set for sequence {genes_list[index][1:12]}: \n"
                    f"forward primer: {f}  Tm: {f_Tm}, GC_content: {f_GC} \n"
                    f"reverse primer: {r}  Tm: {r_Tm}, GC_content: {r_GC}\n"
                    f"Tm difference: {abs(f_Tm - r_Tm)}\n"
                    f"________________________________________________________________________________________\n")
            else:

                print(
                    f"The script found these primers for sequence [ {genes_list[index][1:12]} ] but Tm of the candidate primers was not "
                    f"in the desired range (55-62°C) or diverged more than 4°C in Tm from each other\n "
                    f"forward primer: {f}  Tm: {f_Tm}, GC_content: {f_GC} \n"
                    f"reverse primer: {r}  Tm: {r_Tm}, GC_content: {r_GC} \n"
                    f"The difference is {abs(f_Tm - r_Tm)}°C \n")
                try:
                    f_alter = "NA"
                    r_alter = "NA"

                    # Applying one artificial mismatch to optimize the Tm
                    if f_Tm < 55:
                        f_alter = replacenth(f, "A", "G", 3)  # replaces the 3rd A with G

                    if f_Tm > 62:
                        f_alter = replacenth(f, "G", "A", 3)  # replaces the 3rd G with A

                    if r_Tm < 55:
                        r_alter = replacenth(r, "T", "C", 3)  # replaces the 3rd T with C

                    if r_Tm > 62:
                        r_alter = replacenth(r, "C", "T", 3)  # replaces the 3rd C with T

                    f_alter_Tm, f_alter_GC = tm_gc(f_alter)
                    r_alter_Tm, r_alter_GC = tm_gc(r_alter)

                    if f_alter_GC != 0:
                        print(f"Suggested alternative for forward primer: {f_alter} Tm: {f_alter_Tm} GC: {f_alter_GC}\n")
                    if r_alter_GC != 0:
                        print(f"Suggested alternative for reverse primer: {r_alter} Tm: {r_alter_Tm} GC: {r_alter_GC}\n")
                except:
                    pass

"""
Authors:
	Daniel Prib      	<bytes@mailbox.org>
	Maximilian Leo Huber    <huber.maximilian.leo@gmail.com>
	Saint Fischer 	 	<...>
Version: 1.0
Python Version: 3.11.3
Dependencies: biopython (1.81), numpy (1.25.0)
License: MIT License
"""

import argparse
import os
import time
import csv
import math
from datetime import date
from Bio import SeqIO
from Bio import Seq
import csv
import pathlib
import tempfile
import pandas as pd
import numpy as np
import xml.etree.cElementTree as et
from tabulate import tabulate

# get starting time of program
program_start_time = time.time()

# ----
# util
# ----

def write_csv(output_file, rows):
	with open(output_file, "w") as fp:
		csvwriter = csv.writer(fp, dialect='excel', delimiter=',')
		for row in rows:
			csvwriter.writerow(row)


# ------------
# Saint's Teil
# ------------



# read in arguments, not with sys.argv but with the fancy python thingy
# input: .csv file with all triplets
# flags:
# -o [--output]              set path where output file should be created
# -n [--name]                name output file. Standard value could be time of creation + something else or just time of creation
# -l [--length]              determines the length of the entire area  (I.E: -l 99 => 43bp + Triplet + 43bp). A default value should be set
# -p [--print]               prints out all scored triplets
# -x [--print-only]          Do not create an output file, simply print the results
# -h [--help]                Prints out all flags

# TODO (low priority) check if file has correct extension
# TODO (low priority) implement rest of the flags

print("Parsing arguments...")
parser = argparse.ArgumentParser(description='Input file and flags')
parser.add_argument('-i', required=True, dest='file_path', metavar="FILE", help='file path of the input file (must be a .csv)')
parser.add_argument('-rna', required=True, dest='mrna_path', metavar="FILE", help='mRNA sequence to be scored')
parser.add_argument('-o', required=True, dest='output_path', metavar="DIR", help='folder of the output file. Will be created if it doesn\'t exist yet')
#parser.add_argument('-n', '--name', action='store_const', dest='output_name', help='name of the output file (without extension)', nargs=1)
args = parser.parse_args()

data = SeqIO.read(args.file_path, "fasta")
file = data.id + "_triplets.csv" #output

target = data.seq
target = target.transcribe() #if dna is given: dna to rna
w_seq = target #working seq we will cut this

triplets = ["CCA", "GCA", "UCA", "CAA"]
lfind = -3                                          #Found Index will be added to lfind and 3. for the first loop it dot need the +3 so we offset with -3. This is bc we cut the working sequenz and the cut length need to be readded to new found index

towrite = [["INDEX", "LOCATION", "TRIPLET"],]

def findlocs(w_seq):
    loc = []                                        #Looking for 5UTR, ORF, 3UTR
    loc.append(w_seq.find("AUG"))                   #Find startcodon
    a = w_seq.find("UAA")                           #Find stoppcodon
    if a < w_seq.find("UGA"):
        a = w_seq.find("UGA")
    if a < w_seq.find("UAG"):
        a = w_seq.find("UAG")
    loc.append(a)
    return(loc)

loc = findlocs(target)

def searchhits(w_seq, triplet, lfind):             #Search input for first triplet and append index to list. cut already searched part and repeat
    #print(len(w_seq))                              #debug stuff
    #print(w_seq)                                   #debug stuff
    find = w_seq.find(triplet)
    if find != -1:
        lfind = lfind + 3 + find                     #Find out if found index is before the first start codon or before the first stop codon
        if lfind < loc[0]:
            towrite.append([lfind, "5UTR", triplet])
        else:
            if lfind < loc[1]:
                towrite.append([lfind, "ORF", triplet])
            else:
                towrite.append([lfind, "3UTR", triplet])


        #print(w_seq[find:(find+3)])                #debug stuff

        w_seq = w_seq[find+3::]                     #Cut searched part from the list. Biopython.find only returns index of first match in Sequence ig. Thats why we cut and add found index + lfind
        searchhits(w_seq, triplet, lfind)          #Recursion


for trip in triplets:                              #For each definded triplet we can run searchits
    searchhits(w_seq, trip, lfind)
    w_seq = target
    lfind = -3

write_csv(file, towrite)


# ---------
# Max' Teil
# ---------

if not os.path.exists(file):
    parser.error("The file does not exist!")

today = date.today()

current_path = str(pathlib.Path(__file__).parent.resolve()) + '/'
outpath = current_path + args.output_path + str(today)

if not os.path.exists(outpath):
     os.makedirs(outpath)
    
area_length = 100
half_length = math.ceil(area_length / 2) + 1

difference_scores = []

blastx_exe = str(current_path) + "/ncbi-blast-2.14.0+/bin/blastx"
human_mrna = str(current_path) + "GRCh38_latest_rna.fna"

csv_data = []

print("Blasting...")
with open(str(file), 'r') as csv_file:
    with open(str(args.mrna_path), 'r') as mrna_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        file_counter = 0
        mrna = mrna_file.read()

        if mrna[0] == '>':
            mrna = mrna[mrna.find('\n'):]
            print("removed header line")

        for row in csv_reader:
            csv_data.append(row)
            if row[0] == "INDEX":
                continue
            
            output_file_name = outpath + "/br_" + str(file_counter) + ".xml"
            
            triplet_area = ""
            # TODO (max): test later
            #if len(row[0]) > half_length and len(row[0]) < len(mrna) - half_length:
            
            left_index = int(row[0]) - half_length
            right_index = int(row[0]) + half_length
            triplet_area = mrna[left_index:right_index]
            if (len(triplet_area) < area_length):
                continue
                
            #print(triplet_area)

            tmpfile = open("tmpfileblast.fasta", "w")
            tmpfile.write(str(triplet_area))
            tmpfile.close()
            
            command = "ncbi-blast-2.14.0+/bin/blastn -query tmpfileblast.fasta -subject " + human_mrna + " -out " + output_file_name + " -outfmt 5"
            print(command)
            os.system(command)

            xml_file = et.parse(output_file_name)

            root = xml_file.getroot()

            for diff in root.iter('Hsp_midline'):
                copycount = diff.text.count('|')
                numOfMismatches = len(diff.text) - copycount
                
                difference_scores.append((int(row[0]), numOfMismatches))
            
            file_counter = file_counter + 1

os.system("rm tmpfileblast.fasta")

score_table = [['TRIPLET', 'MISMATCH', "AREA"]]

counter = 1
prev_score = difference_scores[0][0]
for score in difference_scores:
    if not prev_score == score[0]:
        counter = counter + 1
        prev_score = score[0]
    
    score_table.append([score[0], score[1], csv_data[counter][1]])

write_csv("scored.csv", score_table)

# -------------
# Daniel's Teil
# -------------


# 


# diagnostics
program_end_time = time.time()
final_time = program_end_time - program_start_time
print("Program finished! It took " + str(final_time) + " seconds to finish")
# Author: Maximilian Leo Huber
# Date:   20.06.2023

import argparse
import os
import time
import csv
import math
from datetime import date
import pathlib
import tempfile
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO

# get starting time of program
program_start_time = time.time()

current_path = str(pathlib.Path(__file__).parent.resolve()) + '/'

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

check_file_path = current_path + args.file_path
print(check_file_path)
if not os.path.exists(check_file_path):
    parser.error("The file does not exist!")

if not os.path.exists(current_path + args.output_path):
     os.makedirs(current_path + args.output_path)
    
area_length = 100
half_length = math.ceil(area_length / 2) + 1

triplet_score = []
location_score = []

blastx_exe = str(current_path) + "/ncbi-blast-2.14.0+/bin/blastx"
human_mrna = str(current_path) + "GRCh38_latest_rna.fna"
today = date.today()

print("Blasting...")
with open(str(args.file_path), 'r') as csv_file:
    with open(str(args.mrna_path), 'r') as mrna_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        file_counter = 0
        mrna = mrna_file.read()

        if mrna[0] == '>':
            mrna = mrna[mrna.find('\n'):]
            print("removed header line")

        for row in csv_reader:
            if row[0] == "INDEX":
                continue
            
            output_file_name = str(current_path) + str(args.output_path) + "br_" + str(today) + "_" + str(file_counter) + ".xml"

            triplet_area = ""
            # TODO (max): test later
            #if len(row[0]) > half_length and len(row[0]) < len(mrna) - half_length:
            
            left_index = int(row[0]) - half_length
            right_index = int(row[0]) + half_length
            triplet_area = mrna[left_index:right_index]
            if (len(triplet_area) < area_length):
                print("not long enough")
                continue
                
            print(triplet_area)

            query_file = tempfile.NamedTemporaryFile(delete=False)
            query_file.write(triplet_area.encode())
            query_file.close()
            
            print("Blasting line: " + str(file_counter))
            blastx_cline = NcbiblastnCommandline(cmd=blastx_exe,query=query_file.name,subject=human_mrna,out=output_file_name,outfmt=5)
            blastx_cline()

            file_counter = file_counter + 1

# TODO: (high priority) blast triplet areas with all mRNA sequences found in humans
# --- triplet scores ---
#  1. CCA
#  2. GCA
#  3. UCA
#  4. CAA
# --- subobtimal triplets ---
# (5. CUA)
# (6. ACA)
# (7. CCU/CGA)
# (8. CCC)


# TODO (medium priority): score all areas
# Triplet: CCA (I think) > other triplet > ... etc, look at the thing Daniel showed
# Uniqueness: each mismatch = +1
# Region: 3'UTR is better than CDS for example

# TODO: (medium priority) Create .csv file, with columns: Position, Triplet, Uniqueness, Region
# TODO: (medium priority) print results

# diagnostics
program_end_time = time.time()
final_time = program_end_time - program_start_time
print("Program finished! It took " + str(final_time) + " seconds to finish")

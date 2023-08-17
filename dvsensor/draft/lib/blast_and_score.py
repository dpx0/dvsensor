import os
import csv
import xml.etree.cElementTree as et

# blast and score
def blast_and_score(csv_file, mrna_file):  
    csv_file_opened = open(csv_file, "r+")
    csv = csv_file_opened.readlines()
    csv_file_opened.close()
    
    file_counter = 0
    first_row = True

    if not os.path.isdir("blastresults/"):
        os.system("mkdir blastresults/")

    csv_copy = csv
    csv_copy[0] = csv_copy[0][:-1] + ",MISMATCHES,BLASTHITS\n"

    # blast all triplets
    for Index in range(len(csv)):
        if first_row:
            first_row = False
            continue
        
        cells = csv[Index].split(',')

        tmpfile = open("tmpfileblast.fasta", "w")
        tmpfile.write(cells[6])
        tmpfile.close()
        
        output_file_name = "blastresults/br_" + str(file_counter) + ".xml"
        command = "ncbi-blast-2.14.0+/bin/blastn -query tmpfileblast.fasta -subject " + mrna_file + " -out " + output_file_name + " -outfmt 5"
        os.system(command)
        file_counter = file_counter + 1

        xml_file = et.parse(output_file_name)

        root = xml_file.getroot()

        # count mismatches of all triplet areas
        blast_hits = 0
        smallest_diff = 2147483647
        for diff in root.iter('Hsp_midline'):
            blast_hits = blast_hits + 1
            
            if smallest_diff == 0:
                continue 
            
            copycount = diff.text.count('|')
            numOfMismatches = len(diff.text) - copycount
            
            # if an alignment with more mismatches has been found, set it as the biggest difference
            if smallest_diff > numOfMismatches:
                smallest_diff = numOfMismatches
        
        csv_copy[Index] = csv[Index][:-1] + "," + str(smallest_diff) + "," + str(blast_hits) + "\n"
        
    csv_file_opened = open(csv_file, "w")
    for row in csv_copy:
    	csv_file_opened.write(row)
    csv_file_opened.close()

    os.system("rm tmpfileblast.fasta")

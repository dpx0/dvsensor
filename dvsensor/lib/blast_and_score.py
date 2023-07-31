import os
import csv

# blast and score
def blast_and_score(csv_file, full_seq, mrna, area_length, outpath):
    mrna = mrna_file.read()
    file_counter = 0

    # remove header line
    if mrna[0] == '>':        
        mrna = mrna[mrna.find('\n'):]
        print("removed header line")
    
    # blast all triplets
    for row in csv_file:
        # skip column names
        if row[0] == "INDEX":
            continue
        
        # get triplet area
        half_length = area_length / 2   
         
        left_index = int(row[0]) - half_length
        right_index = int(row[0]) + half_length
        
        # if any triplet area goes out of bounds, skip
        if left_index < 0 or right_index >= len(full_seq):
        	continue
        
        triplet_area = full_seq[left_index:right_index]

        tmpfile = open("tmpfileblast.fasta", "w")
        tmpfile.write(str(triplet_area))
        tmpfile.close()
            
        output_file_name = outpath + "/br_" + str(file_counter) + ".xml"
        command = "ncbi-blast-2.14.0+/bin/blastn -query tmpfileblast.fasta -subject " + mrna + " -out " + output_file_name + " -outfmt 5"
        os.system(command)
        file_counter = file_counter + 1

        xml_file = et.parse(output_file_name)

        root = xml_file.getroot()

        # count mismatches of all triplet areas
        biggest_diff = -1
        for diff in root.iter('Hsp_midline'):
            copycount = diff.text.count('|')
            numOfMismatches = len(diff.text) - copycount
            
            # if a perfect alignment has been found, skip this triplet area
            if numOfMismatches == 0:
            	biggest_diff = -1
            	break
            
            # if an alignment with more mismatches has been found, set it as the biggest difference
            if biggest_diff < numOfMismatches:
                biggest_diff = numOfMismatches
        
        # Skip triplet areas which have perfect alignments
        if biggest_diff == -1:
        	continue
        
        csv_file.append(current_score)

    os.system("rm tmpfileblast.fasta")
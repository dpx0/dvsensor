from Bio import SeqIO
from Bio import Seq
import csv
import argparse

parser = argparse.ArgumentParser(description='Process .FASTA single DNA/RNA into possible DART VADER targets') #übergabe der Argumente
parser.add_argument('filename')
args = parser.parse_args()

inp = args.filename
data = SeqIO.read(inp, "fasta")
file = data.id + ".csv" #output

target = data.seq
target = target.transcribe() #if dna is given: dna to rna
w_seq = target #working seq we will cut this

tripplets = ["CCA", "GCA", "UCA", "CAA"]
lfind = -3                                          #Found Index will be added to lfind and 3. for the first loop it dot need the +3 so we offset with -3. This is bc we cut the working sequenz and the cut length need to be readded to new found index



towrite = [["INDEX", "LOCATION", "TRIPPLET"],]

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

def searchhits(w_seq, tripplet, lfind):             #Search input for first tripplet and append index to list. cut already searched part and repeat
    #print(len(w_seq))                              #debug stuff
    #print(w_seq)                                   #debug stuff
    find = w_seq.find(tripplet)
    if find != -1:
        lfind = lfind + 3 + find                     #Find out if found index is before the first start codon or before the first stop codon
        if lfind < loc[0]:
            towrite.append([lfind, "5UTR", tripplet])
        else:
            if lfind < loc[1]:
                towrite.append([lfind, "ORF", tripplet])
            else:
                towrite.append([lfind, "3UTR", tripplet])


        #print(w_seq[find:(find+3)])                #debug stuff

        w_seq = w_seq[find+3::]                     #Cut searched part from the list. Biopython.find only returns index of first match in Sequence ig. Thats why we cut and add found index + lfind
        searchhits(w_seq, tripplet, lfind)          #Recursion




for trip in tripplets:                              #For each definded tripplet we can run searchits
    searchhits(w_seq, trip, lfind)
    w_seq = target
    lfind = -3


with open(file, 'w', newline="\n") as f:            #Save as .csv
    mywriter = csv.writer(f, dialect='excel', delimiter=',')
    mywriter.writerows(towrite)

#print(loc)                                         #debug stuff


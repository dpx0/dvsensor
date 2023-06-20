from Bio import SeqIO
from Bio import Seq
import csv
import numpy as np

#Diese Datei sucht die angegebenen Tripplets und gibt die stellen zurück. EIngabe soll eine dna/rna .fasta und die gesuchten tripplets sein


inp = "IDH1_ref.fasta" #input soll später mit arg übergeben werden

data = SeqIO.read(inp, "fasta")
target = data.seq
target = target.transcribe() #dna zu rna
w_seq = target #working seq. I will cut this
hit = []
tripplet = "CCA"
tripplets = ["CCA", "GCA", "UCA", "CAA"]
lfind = -3 #jeder find wird mit lfind und 3 addiert um die stelle der seq zu finden
file = data.id + ".csv"
loc = []

towrite = [["INDEX", "LOCATION", "TRIPPLET"],]

def findlocs(w_seq): #Funktion sucht nach 5UTR, ORF, 3UTR
    loc.append(w_seq.find("AUG")) #Finde das startcodon
    a = w_seq.find("UAA") #Finde das stoppcodon
    if a < w_seq.find("UGA"):
        a = w_seq.find("UGA")
    if a < w_seq.find("UAG"):
        a = w_seq.find("UAG")
    loc.append(a)
    return(loc)

loc = findlocs(target)

def searchhits(w_seq, tripplet, lfind): #Funktion sucht eingabesequenz nach der erstgelegenen trippletstelle, gibt sie an hit ab und schneidet die seq. wiederholung bis keine seq mehr da
    #print(len(w_seq))
    print(w_seq)
    find = w_seq.find(tripplet)
    if find != -1:
        lfind = lfind + 3 + find
        if lfind < loc[0]:
            towrite.append([lfind, "5UTR", tripplet])
        else:
            if lfind < loc[1]:
                towrite.append([lfind, "ORF", tripplet])
            else:
                towrite.append([lfind, "3UTR", tripplet])
        #towrite.append([lfind, loc, tripplet])
        hit.append((lfind))
        print(hit)
        print(w_seq[find:(find+3)])
        w_seq = w_seq[find+3::] #Die liste wird stückweise kleiner find weiß wo die nächste seq ist
        #print(len(w_seq)+hit)
        searchhits(w_seq, tripplet, lfind)




for trip in tripplets: #Jedes Angegebene Tripplet wird einmal durchlaufen
    searchhits(w_seq, trip, lfind)
    w_seq = target
    lfind = -3

for index in range(len(hit)): #überprüft ob die Ergebnisse stimmen, ausgabe aller stellen, die hit gefunden hat
    print(target[hit[index]:hit[index]+3])


with open('sample.csv', 'w', newline="\n") as f: #das ganze als csv speichern
    mywriter = csv.writer(f, dialect='excel', delimiter=',')
    mywriter.writerows(towrite)

print(loc)

#np.savetxt('sample.csv',arr, fmt = '%d', delimiter=",")

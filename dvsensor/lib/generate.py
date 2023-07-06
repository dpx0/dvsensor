"""
Authors:
	Daniel Prib      		<bytes@mailbox.org>
	Maximilian Leo Huber    <...>
	Saint Fischer 	 		<...>
Version: 1.0
Python Version: 3.11.3
Dependencies: biopython (1.81), numpy (1.25.0)
License: MIT License
"""

import re
from Bio.Seq import Seq


MS2_HAIRPIN = "ACAUGAGGAUCACCCAUGU"


def generate_sensor_sequences(sequence, triplets, ntleft, ntright):
	for entry in triplets:
		start = entry[1] - ntleft
		stop = entry[1] + 3 + ntright
		upstream = str(Seq(sequence[start: entry[1]]).reverse_complement())
		downstream = str(Seq(sequence[entry[1] + 3: stop]).reverse_complement())
		# the order of upstream and downstream is reversed, because these are reverse complements
		sensorSeq = Seq(process_part(downstream) + "UAG" + process_part(upstream)).back_transcribe()
		triggerSeq = Seq(sequence[start : stop]).back_transcribe()[::-1]
		entry.extend([str(start) + "-" + str(stop), str(sensorSeq), str(triggerSeq)])

	return triplets


def process_part(seq):

	# insert MS2-hairpins
	seq = seq[::-1]
	seq = seq[:24] + MS2_HAIRPIN[::-1] + seq[31:]

	## edit in-frame stopcodons (reversed)
	for codon in ('GAU', 'AAU', 'AGU'):
		positions = [match.start() for match in re.finditer(codon, seq)]
		positions = [pos for pos in positions if (pos % 3) == 0]

		for pos in positions:
			seq = seq[:pos] + {'GAU': 'GAC',
							   'AAU': 'ACU',
							   'AGU': 'AGG'}.get(codon) + seq[(pos + 3):]

	return seq[::-1]


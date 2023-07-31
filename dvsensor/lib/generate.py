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
nStopEdits = 0


def generate_sensor_sequences(sequence, triplets, ntleft, ntright):
	global nStopEdits

	for entry in triplets:
		nStopEdits = 0
		start = entry[1] - ntleft
		stop = entry[1] + 3 + ntright

		upstreamRC = str(Seq(sequence[start: entry[1]]).reverse_complement())
		downstreamRC = str(Seq(sequence[entry[1] + 3: stop]).reverse_complement())

		# the order of upstream and downstream is reversed, because these are reverse complements
		sensorSeq = Seq(process_part(downstreamRC, 'LEFT') +
						"UAG" +
						process_part(upstreamRC, 'RIGHT')).back_transcribe()

		triggerSeq = Seq(sequence[start:stop]).back_transcribe()[::-1]  # 3'-> 5'
		entry.extend([str(start) + "-" + str(stop), str(sensorSeq), str(triggerSeq),
					  str(triggerSeq[::-1]), nStopEdits])

	return triplets


def process_part(seq, side):
	global nStopEdits

	# insert MS2-hairpins
	seq = seq[::-1]
	if side == 'LEFT':
		seq = seq[:24] + MS2_HAIRPIN[::-1] + seq[31:]
	elif side == 'RIGHT':
		seq = seq[:18] + MS2_HAIRPIN[::-1] + seq[25:]

	## edit in-frame stopcodons (reversed)
	for codon in ('GAU', 'AAU', 'AGU'):
		positions = [match.start() for match in re.finditer(codon, seq)]
		positions = [pos for pos in positions if (pos % 3) == 0]
		nStopEdits += len(positions)

		for pos in positions:
			seq = seq[:pos] + {'GAU': 'GAC',
							   'AAU': 'ACU',
							   'AGU': 'AGG'}.get(codon) + seq[(pos + 3):]

	return seq[::-1]


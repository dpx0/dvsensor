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

from Bio import SeqIO
from csv import writer
import re


def locate_triplets(input_file, triplets, regions):

	inputSeq = str(SeqIO.read(input_file, "fasta").seq.transcribe())  # assume cDNA as input sequence
	startPos = inputSeq.find("AUG")

	# find all potential stopcodons
	allStops = []
	for c in ("UAA", "UAG", "UGA"):
		allStops.extend([match.start() for match in re.finditer(c, inputSeq)])

	# filter in-frame stopcodons and choose the first occurence
	stopPos = min([p for p in allStops if (p > startPos and (p - startPos) % 3 == 0)])

	results = {region: [] for region in regions}

	for triplet in triplets:
		positions = [match.start() for match in re.finditer(triplet, inputSeq)]

		for pos in positions:
			if pos < startPos:
				if '5UTR' in regions: results['5UTR'].append((pos, triplet))

			elif pos > stopPos:
				if '3UTR' in regions: results['3UTR'].append((pos, triplet))

			elif 'CDS' in regions: results['CDS'].append((pos, triplet))

	return results


def write_output(output_file, results):

	pass




import logging
import re
import io
import subprocess
import Bio.Blast.NCBIXML
from collections import namedtuple
from Bio.Seq import Seq

TranscriptRegions = namedtuple('TranscriptRegions', ['coord_5UTR', 'coord_CDS', 'coord_3UTR'])
Triplet = namedtuple('Triplet', ['region', 'pos', 'kind'])
SensorData = namedtuple('SensorData', ['sensor_seq', 'trigger_seq', 'num_edits'])
BlastHit = namedtuple('BlastHit', ['accession', 'expect', 'score', 'q_coverage'])


def find_transcript_regions(rna_sequence: str) -> TranscriptRegions | None:
	start_pos = rna_sequence.find('AUG')
	if start_pos == -1:  # 'AUG' not found
		return None
	# find all potential stop-codons
	all_stops = []
	for c in ("UAA", "UAG", "UGA"):
		all_stops.extend([match.start() for match in re.finditer(c, rna_sequence)])
	# filter for downstream and in-frame stop-codons and choose the first occurrence
	try:
		stop_pos = min([p for p in all_stops if (p > start_pos and (p - start_pos) % 3 == 0)])
	except ValueError:  # no valid stop-codon found
		return None

	return TranscriptRegions(coord_5UTR=(0, start_pos),
							 coord_CDS=(start_pos, stop_pos+3),
							 coord_3UTR=(stop_pos+3, len(rna_sequence)))


def find_triplets(sequence: str, target_triplets: tuple[str], target_regions: tuple[str],
				  coord_CDS: tuple[int, int]) -> list[Triplet]:
	results = []
	for triplet in target_triplets:
		positions = [match.start() for match in re.finditer(triplet, sequence)]
		for pos in positions:
			if pos < coord_CDS[0] and '5UTR' in target_regions:
				results.append(Triplet(region='5UTR',
									   pos=pos,
									   kind=triplet))
			elif pos > coord_CDS[1] and '3UTR' in target_regions:
				results.append(Triplet(region='3UTR',
									   pos=pos,
									   kind=triplet))
			elif 'CDS' in target_regions:
				results.append(Triplet(region='CDS',
									   pos=pos,
									   kind=triplet))
	return results


def generate_sensor(sequence: str, triplet_pos: int) -> SensorData | None:
	ntleft = 48
	ntright = 48
	sensor_start = triplet_pos - ntleft
	sensor_end = triplet_pos + 3 + ntright
	# dismiss sensors with only partial lengths (<123 bp total)
	if sensor_start < 0 or sensor_end > len(sequence) - 1:
		return

	trigger_seq = sequence[sensor_start: sensor_end]  # 5' -> 3'
	sensor_seq = str(Seq(trigger_seq).reverse_complement_rna())  # 5' -> 3'

	# insert central stop codon + MS2 hairpins
	sensor_seq = insert_hairpins(f'{sensor_seq[:ntleft]}UAG{sensor_seq[(ntleft + 3):]}', ntleft)

	# edit in-frame start- and stop-codons
	left_part, left_edits = edit_in_frame_codons(sensor_seq[:(ntleft + 12)])
	right_part, right_edits = edit_in_frame_codons(sensor_seq[(ntleft + 15):])
	sensor_seq = f'{left_part}UAG{right_part}'

	return SensorData(sensor_seq=sensor_seq,
					  trigger_seq=trigger_seq,
					  num_edits=(left_edits + right_edits))


def edit_in_frame_codons(sequence: str) -> tuple[str, int]:
	num_edits = 0
	substitutions = {'UAG': 'CAG', 'UGA': 'GGA',
					 'UAA': 'UCA', 'AUG': 'AUU'}

	for codon in ('UAG', 'UGA', 'UAA', 'AUG'):
		matches = [match.start() for match in re.finditer(codon, sequence)]
		in_frame_matches = [match for match in matches if (match % 3) == 0]
		num_edits += len(in_frame_matches)

		for match in in_frame_matches:
			sequence = sequence[:match] + substitutions.get(codon) + sequence[(match + 3):]
	return sequence, num_edits


def insert_hairpins(sequence: str, ntleft: int) -> str:
	MS2_HAIRPIN = 'ACAUGAGGAUCACCCAUGU'  # MS2 bacteriophage hairpin sequence
	hairpin_lpos = ntleft - 31
	hairpin_rpos = ntleft + 26
	return (
		f'{sequence[:hairpin_lpos]}{MS2_HAIRPIN}{sequence[(hairpin_lpos + 7):ntleft]}'
		f'{sequence[ntleft:(ntleft + 3)]}'
		f'{sequence[(ntleft + 3):hairpin_rpos]}{MS2_HAIRPIN}{sequence[(hairpin_rpos + 7):]}'
	)


def blast_analysis(query_sequence: str,
				   query_accession: str,
				   only_overlapping: bool,
				   blast_variant: str,
				   db_name: str,
				   word_size: str,
				   evalue: str,
				   taxids: str,
				   perc_identity: str,
				   qcov_hsp_perc: str,
				   db_path: str) -> list[BlastHit]:
	query_accession = query_accession.split('.')[0]
	blastcmd = ["blastn",
				"-task", blast_variant,
				"-db", db_name,
				"-word_size", word_size,
				"-evalue", evalue,
				"-outfmt", "5"]
	if taxids:
		blastcmd.extend(["-taxids", taxids])
	if perc_identity:
		blastcmd.extend(["-perc_identity", perc_identity])
	if qcov_hsp_perc:
		blastcmd.extend(["-qcov_hsp_perc", qcov_hsp_perc])

	blast_result = subprocess.run(blastcmd, input=query_sequence,
								  capture_output=True, text=True,
								  env={'BLASTDB': db_path})
	try:
		blast_record = Bio.Blast.NCBIXML.read(io.StringIO(blast_result.stdout))
	except ValueError:
		logging.warning("empty blast record")
		return []

	blast_hits = []
	query_len = len(query_sequence)
	ntleft = (query_len - 3)/2
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if not alignment.accession.split('.')[0] == query_accession:
				if only_overlapping and (hsp.query_start > (ntleft + 1) or hsp.query_end < (ntleft + 3)):
					continue
				blast_hits.append(BlastHit(accession=alignment.accession,
										   expect=hsp.expect,
										   score=hsp.score,
										   q_coverage=(hsp.query_end - hsp.query_start) / query_len))
	return blast_hits

import logging
import re
import io
import threading
import subprocess
from typing import Any, Callable
import Bio.Blast.NCBIXML
from interface import SignalHandler, Signal
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def generate_sensors_job(job_data: dict[str, Any],
						 stop_event: threading.Event,
						 signals_handler: SignalHandler) -> None:

	sequence_data = job_data['sequence_data']
	options = job_data['options']
	sequence: str = sequence_data['sequence']

	signals_handler.register_signal('ideogram_data_ready')
	transcript_regions: dict[str, tuple[int, int]] = find_transcript_regions(sequence)
	ideogram_data: dict[str, int | float] = calculate_ideogram_data(transcript_regions)
	signals_handler.send_signal('ideogram_data_ready', args=ideogram_data)

	ui_ready: Signal = signals_handler.request_signal('ui_ready')
	ui_functions: dict[str, Callable] = ui_ready.wait()
	update_progress: Callable = ui_functions['update_progress']
	set_status_finished: Callable = ui_functions['set_status_finished']
	add_rows: Callable = ui_functions['add_rows']

	CDS_start, CDS_end = transcript_regions['CDS']
	triplets: list[tuple[str, int, str]] = \
		find_triplets(sequence,
					  include_triplets=options['triplets'],
					  include_regions=options['regions'],
					  CDS_start=CDS_start,
					  CDS_end=CDS_end)

	progress_step: float = (1.0 / len(triplets) if triplets else 100.0)
	_finished = False
	while stop_event and not stop_event.is_set():
		if not triplets:
			_finished = True
			break

		region, position, triplet = triplets.pop()
		result: tuple[str, str, int] = generate_sensor(sequence, position)

		if not result:
			continue
		sensor, trigger, edits = result
		blast_off_targets: list[str] = run_blast_analysis(trigger, job_options)
		print(";".join(blast_off_targets))
		sensor_entry = {
			'position': position + 1,  # nucleotide position counting from 1
			'triplet': triplet,
			'region': region,
			'range': f'{(position + 1) - 48}-{(position + 1) + 50}',
			'percent_gc': round(GC(sensor), 1),
			'n_edits': edits,
			'off_targets': ";".join(blast_off_targets),
			'sensor': sensor,
			'trigger': trigger
		}
		add_rows([sensor_entry])
		update_progress(progress_step)

	if _finished:
		set_status_finished()
		logging.debug('job finished')
	else:
		logging.debug('job did not finish')


def find_transcript_regions(sequence: str) -> dict[str, tuple[int, int]]:
	start_pos = sequence.find('AUG')
	# find all potential stopcodons
	all_stops = []
	for c in ("UAA", "UAG", "UGA"):
		all_stops.extend([match.start() for match in re.finditer(c, sequence)])
	# filter for in-frame stopcodons and choose the first occurence
	stop_pos = min([p for p in all_stops if (p > start_pos and (p - start_pos) % 3 == 0)])
	return {'5UTR': (0, start_pos),
			'CDS': (start_pos, stop_pos+3),
			'3UTR': (stop_pos+3, len(sequence))}


def calculate_ideogram_data(transcript_regions: dict[str, tuple[int, int]]) -> dict[str, int | float]:
	len_5UTR = transcript_regions['5UTR'][1] - transcript_regions['5UTR'][0]
	len_CDS = transcript_regions['CDS'][1] - transcript_regions['CDS'][0]
	len_3UTR = transcript_regions['3UTR'][1] - transcript_regions['3UTR'][0]
	len_tot = len_5UTR + len_CDS + len_3UTR
	percent_5UTR = len_5UTR / len_tot * 100
	percent_CDS = len_CDS / len_tot * 100
	percent_3UTR = 100.0 - (percent_5UTR + percent_CDS)

	return {
		'len_5UTR': len_5UTR,
		'len_CDS': len_CDS,
		'len_3UTR': len_3UTR,
		'percent_5UTR': percent_5UTR,
		'percent_CDS': percent_CDS,
		'percent_3UTR': percent_3UTR
	}


def find_triplets(sequence: str, include_triplets: tuple, include_regions: tuple,
				  CDS_start: int, CDS_end: int) -> list[tuple[str, int, str]]:
	results = []
	for triplet in include_triplets:
		positions = [match.start() for match in re.finditer(triplet, sequence)]
		for pos in positions:
			if pos < CDS_start:
				if '5UTR' in include_regions:
					results.append(('5UTR', pos, triplet))
			elif pos > CDS_end:
				if '3UTR' in include_regions:
					results.append(('3UTR', pos, triplet))
			elif 'CDS' in include_regions:
				results.append(('CDS', pos, triplet))
	return results


def generate_sensor(sequence: str, triplet_position: int) -> tuple[str, str, int] | None:

	def _edit_codons(sequence_part: str) -> tuple[str, int]:
		edits = 0
		for codon in ('UAG', 'UGA', 'UAA', 'AUG'):
			matches = [match.start() for match in re.finditer(codon, sequence_part)]
			in_frame_matches = [match for match in matches if (match % 3) == 0]
			edits += len(in_frame_matches)

			for match in in_frame_matches:
				sequence_part = sequence_part[:match] + \
								{'UAG': 'CAG',
								 'UGA': 'GGA',
								 'UAA': 'UCA',
								 'AUG': 'AUU'}.get(codon) + sequence_part[(match + 3):]
		return sequence_part, edits

	ntleft = 48
	ntright = 48
	sensor_start = triplet_position - ntleft
	sensor_end = triplet_position + 3 + ntright
	# dismiss sensors with only partial lengths (<123 bp total)
	if sensor_start < 0 or sensor_end > len(sequence) - 1:
		return

	trigger_region = sequence[sensor_start: sensor_end]  # 5' -> 3'
	sensor = str(Seq(trigger_region).reverse_complement_rna())  # 5' -> 3'
	# insert central stop codon
	sensor = f'{sensor[:ntleft]}UAG{sensor[(ntleft + 3):]}'

	# insert MS2 hairpins
	MS2_HAIRPIN = 'ACAUGAGGAUCACCCAUGU'
	left_MS2_start = ntleft - 31
	right_MS2_start = ntleft + 26
	sensor = f'{sensor[:left_MS2_start]}{MS2_HAIRPIN}{sensor[(left_MS2_start + 7):ntleft]}' + \
			 f'{sensor[ntleft:(ntleft + 3)]}' + \
			 f'{sensor[(ntleft + 3):right_MS2_start]}{MS2_HAIRPIN}{sensor[(right_MS2_start + 7):]}'

	# edit in-frame start- and stop-codons
	left_part, left_edits = _edit_codons(sensor[:(ntleft + 12)])
	right_part, right_edits = _edit_codons(sensor[(ntleft + 15):])
	sensor = f'{left_part}UAG{right_part}'

	return sensor, trigger_region, (left_edits + right_edits)


def run_blast_analysis(query_sequence: str, **kwargs) -> list[str]:
	blastcmd = [
		"blastn", "-db", kwargs['db'], "-taxids", kwargs['taxids'],
		"-word_size", kwargs['word_size'], "-outfmt", "5",
		"-evalue", kwargs['e_value']
	]
	query_acc_num = kwargs['accession_num'].split('.')[0]

	blast_result = subprocess.run(blastcmd, capture_output=True, text=True, env={'BLASTDB': kwargs['blastdb']},
								  input=query_sequence)
	blast_record = Bio.Blast.NCBIXML.read(io.StringIO(blast_result.stdout))
	off_targets = []
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			#if hsp.expect <0.01:
			print(hsp.expect)
			if not alignment.accession.split('.')[0] == query_acc_num:
				off_targets.append(f'{alignment.accession}(E:{hsp.expect};Score:{hsp.score})')
				# print('****Alignment****')
				# print('sequence:', alignment.title)
				# print('accession: ', alignment.accession)
				# print('length:', alignment.length)
				# print('score:', hsp.score)
				# print('gaps:', hsp.gaps)
				# print('e value:', hsp.expect)
				# print(hsp.query[0:90] + '...')
				# print(hsp.match[0:90] + '...')
				# print(hsp.sbjct[0:90] + '...')
		print()
	return off_targets

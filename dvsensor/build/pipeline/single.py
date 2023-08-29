from typing import Optional
from Bio.Seq import Seq
import re
import time


def heavy_computation():
	...
	# time.sleep(1)


def analyze_single_sequence(task_handler, task_options: dict) -> None:

	ui_connection = task_handler.ui_connection
	data = task_handler.data
	sequence = data.rna_data.sequence
	# if task_options['file_format'] == 'FASTA':
	transcript_regions = find_transcript_regions(sequence)
	data.update('len_5UTR', transcript_regions['5UTR'][1] - transcript_regions['5UTR'][0])
	data.update('len_CDS', transcript_regions['CDS'][1] - transcript_regions['CDS'][0])
	data.update('len_3UTR', transcript_regions['3UTR'][1] - transcript_regions['3UTR'][0])

	include_triplets = tuple(triplet for triplet in task_options['triplet_settings'].keys() if
				task_options['triplet_settings'][triplet])
	include_regions = tuple(region for region in task_options['regions_settings'].keys() if
				task_options['regions_settings'][region])
	start_pos, stop_pos = transcript_regions['CDS']
	found_triplets = find_triplets(sequence, include_triplets, include_regions, start_pos, stop_pos)
	data.update('num_triplets', len(found_triplets))

	ui_connection.page_render_complete.wait()
	progress_step = 1.0 / len(found_triplets)
	progress_bar = ui_connection.get_element('progress_bar')
	while task_handler.thread:
		if not found_triplets:
			ui_connection.call('set_status_finished')
			task_handler.terminate()
			break

		region, position, triplet = found_triplets.pop()
		result = generate_sensor(sequence, position)
		heavy_computation()  # <---------------------- dummy function
		if not result:
			continue
		sensor, trigger, edits = result
		sensor_entry = {
			'position': position + 1,  # nucleotide position counting from 1
			'triplet': triplet,
			'region': region,
			'range': f'{(position + 1) - 48}-{(position + 1) + 50}',
			'percent_gc': calc_gc_content(sensor),
			'n_edits': edits,
			'off_targets': '-',
			'sensor': sensor,
			'trigger': trigger
		}
		ui_connection.call('add_rows', row_data=[sensor_entry])
		ui_connection.call('update_progress', step=progress_step)
		if round(progress_bar.value * 100) >= 100:
			progress_bar.value = 1.0
	print("thread finished")


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


def find_triplets(sequence: str, include_triplets: tuple, include_regions: tuple,
				  start_pos: int, stop_pos: int) -> list[tuple[str, int, str]]:

	results = []
	for triplet in include_triplets:
		positions = [match.start() for match in re.finditer(triplet, sequence)]

		for pos in positions:
			if pos < start_pos:
				if '5UTR' in include_regions:
					results.append(('5UTR', pos, triplet))

			elif pos > stop_pos:
				if '3UTR' in include_regions:
					results.append(('3UTR', pos, triplet))

			elif 'CDS' in include_regions:
				results.append(('CDS', pos, triplet))
	return results


def generate_sensor(sequence: str, triplet_position: int) -> Optional[tuple[str, str, int]]:

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
	sensor = str(Seq(trigger_region).reverse_complement())  # 5' -> 3'
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


def calc_gc_content(sequence: str) -> float:
	sequence = list(sequence)
	gc_content = (sequence.count('C') + sequence.count('G')) / len(sequence)
	return round(gc_content * 100, 1)

	# num_gc = 0
	# for b in sequence:
	# 	if b == 'C' or b == 'G':
	# 		num_gc += 1
	# return round((num_gc / len(sequence)) * 100, 1)

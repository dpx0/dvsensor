from pprint import pprint
import re
import asyncio


async def analyze_single_sequence(task_handler, task_options) -> None:

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
	pprint(found_triplets)

	# print(sequence[transcript_regions['5UTR'][0] : transcript_regions['5UTR'][1]])
	# print(sequence[transcript_regions['CDS'][0]: transcript_regions['CDS'][1]])
	# print(sequence[transcript_regions['3UTR'][0]: transcript_regions['3UTR'][1]])


	# row_data = [
	# 	{'position': '7127',
	# 	 'triplet': 'CCA',
	# 	 'region': '3UTR',
	# 	 'range': '7079-7178',
	# 	 'percent_gc': '52,8%',
	# 	 'n_stop_codons': '0',
	# 	 'off_targets': '0'},
	# 	{'position': '5788',
	# 	 'triplet': 'CCA',
	# 	 'region': '3UTR',
	# 	 'range': '5740-5839',
	# 	 'percent_gc': '51,2%',
	# 	 'n_stop_codons': '0',
	# 	 'off_targets': '0'},
	# ]
	#
	# progress_step = 1.0 / len(row_data)
	# progress_bar = ui_connection.get_element('progress_bar')
	# while True:
	# 	await asyncio.sleep(2)
	# 	if not row_data:
	# 		break
	# 	ui_connection.call('add_rows', [row_data.pop()])
	#
	# 	ui_connection.call('update_progress', step=progress_step)
	# 	print(progress_bar.value)
	# 	if round(progress_bar.value * 100) >= 100:
	# 		print("FINISHED")
	# 		ui_connection.call('set_status_finished')
	# 		task_handler.terminate()
	# 		break


def find_transcript_regions(sequence: str) -> dict[str, tuple[int, int]]:
	start_pos = sequence.find('AUG')

	# find all potential stopcodons
	all_stops = []
	for c in ("UAA", "UAG", "UGA"):
		all_stops.extend([match.start() for match in re.finditer(c, sequence)])

	# filter for in-frame stopcodons and choose the first occurence
	stop_pos = min([p for p in all_stops if (p > start_pos and (p - start_pos) % 3 == 0)])

	return {'5UTR': (0, start_pos),
			'CDS': (start_pos+1, stop_pos+3),
			'3UTR': (stop_pos+4, len(sequence))}


def find_triplets(sequence: str, include_triplets: tuple, include_regions: tuple,
				  start_pos: int, stop_pos: int) -> dict[str: list[tuple[int, str]]]:

	results = {region: [] for region in include_regions}

	for triplet in include_triplets:
		positions = [match.start() for match in re.finditer(triplet, sequence)]

		for pos in positions:
			if pos < start_pos:
				if '5UTR' in include_regions:
					results['5UTR'].append((pos, triplet))

			elif pos > stop_pos:
				if '3UTR' in include_regions:
					results['3UTR'].append((pos, triplet))

			elif 'CDS' in include_regions:
				results['CDS'].append((pos, triplet))

	return results

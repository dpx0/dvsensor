import threading
import logging
import jobs.pipeline as pl
from typing import Any, Callable
from Bio.SeqUtils import gc_fraction
from interface import SignalHandler, Signal


class JobRuntimeError(Exception):
	def __init__(self, message) -> None:
		self.message = message
		super().__init__(self.message)


def generate_sensors_job(job_data: dict[str, Any],
						 stop_event: threading.Event,
						 signals_handler: SignalHandler) -> None:
	try:
		sequence: str = job_data['sequence_data']['sequence']
		sequence_id: str = job_data['sequence_data']['id']
		target_triplets: tuple[str] = job_data['options']['triplets']
		target_regions: tuple[str] = job_data['options']['regions']
		blast_options: dict | None = job_data['options']['blast']
		if not blast_options['use_blast']:
			blast_options = None
		else:
			del blast_options['use_blast']
	except KeyError as e:
		logging.error(f"Job Error: missing key {e} in job_data")
		return

	required_ui_fn = ['update_progress',
					  'set_status_finished',
					  'add_rows']
	try:
		transcript_regions: pl.TranscriptRegions = pl.find_transcript_regions(sequence)
		if not transcript_regions:
			raise JobRuntimeError('transcript regions not found')
		provide_ideogram_data(calc_ideogram_data(transcript_regions), signals_handler)
		# ui_functions waits for a 'ui_ready' signal and needs to be called after send_ideogram_data
		ui_functions: dict[str, Callable] = get_ui_functions(signals_handler)
		for fn in required_ui_fn:
			if fn not in ui_functions:
				raise JobRuntimeError(f'required ui function {fn} not received')

		triplets: list[pl.Triplet] = pl.find_triplets(sequence=sequence,
													  target_triplets=target_triplets,
													  target_regions=target_regions,
													  coord_CDS=transcript_regions.coord_CDS)
	except JobRuntimeError as e:
		logging.error(f'Job Error: {e.message}')
		return

	# TODO: supply these keys in blast_options
	# ----------------------------------------
	if blast_options:
		blast_options['db_name'] = 'refseq_select_rna'
		blast_options['db_path'] = './test/blastdb'
		blast_options['taxids'] = '9606'
		blast_options['blast_variant'] = 'blastn'
	# ----------------------------------------

	try:
		if mainloop(stop_event=stop_event,
					sequence=sequence,
					sequence_id=sequence_id,
					triplets=triplets,
					blast_options=blast_options,
					add_rows=ui_functions['add_rows'],
					update_progress=ui_functions['update_progress']):
			ui_functions['set_status_finished']()
			logging.debug("job finished")
		else:
			logging.debug("job did not finish")
	except JobRuntimeError as e:
		logging.error(f'Job Runtime Error: {e.message}')


def mainloop(stop_event: threading.Event, sequence: str, sequence_id: str, triplets: list[pl.Triplet],
			 blast_options: dict | None, add_rows: Callable, update_progress: Callable) -> bool:

	if not triplets:
		update_progress(100.0)
		return True

	progress_step: float = 1.0 / len(triplets)
	potential_off_targets: str = 'N.A.'
	while stop_event and not stop_event.is_set():
		if not triplets:
			return True

		triplet = triplets.pop()
		sensor_data: pl.SensorData = pl.generate_sensor(sequence, triplet.pos)
		if not sensor_data:
			continue  # sensors with only partial lengths (<123 bp total) are dismissed

		if blast_options:
			blast_results: list[pl.BlastHit] = pl.blast_analysis(query_sequence=sensor_data.trigger_seq,
																 query_accession=sequence_id,
																 **blast_options)
			potential_off_targets = ";".join([
				f'{hit.accession}(E:{hit.expect};Score:{hit.score};Cov:{hit.q_coverage * 100:.0f}%)'
				for hit in blast_results])

		sensor_entry = {
			'position': triplet.pos + 1,  # first nucleotide position of triplet in transcript, counting from 1
			'triplet': triplet.kind,
			'region': triplet.region,
			'range': f'{(triplet.pos + 1) - 48}-{(triplet.pos + 1) + 50}',
			'percent_gc': f'{gc_fraction(sensor_data.sensor_seq) * 100:.1f}',
			'n_edits': sensor_data.num_edits,
			'pot_off_targets': potential_off_targets,
			'sensor': sensor_data.sensor_seq,
			'trigger': sensor_data.trigger_seq
		}
		add_rows([sensor_entry])
		update_progress(progress_step)
	return False


def get_ui_functions(signals_handler: SignalHandler) -> dict[str, Callable]:
	ui_ready: Signal = signals_handler.request_signal('ui_ready')
	ui_functions: dict[str, Callable] = ui_ready.wait()
	return ui_functions


def calc_ideogram_data(transcript_regions: pl.TranscriptRegions) -> dict[str, int | float]:
	len_5UTR = transcript_regions.coord_5UTR[1] - transcript_regions.coord_5UTR[0]
	len_CDS = transcript_regions.coord_CDS[1] - transcript_regions.coord_CDS[0]
	len_3UTR = transcript_regions.coord_3UTR[1] - transcript_regions.coord_3UTR[0]
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


def provide_ideogram_data(ideogram_data: dict[str, int | float],
						  signals_handler: SignalHandler) -> None:
	signals_handler.register_signal('ideogram_data_ready')
	signals_handler.send_signal('ideogram_data_ready', args=ideogram_data)

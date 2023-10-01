import os
import logging
from pathlib import Path
from nicegui import ui
from typing import Any, Callable
from interface import Controller
from jobs import generate_sensors_job
from .style import Colors
from .pagebuilder import page_builder
from .ui_elements import back_button, show_page_not_available, show_dialog_box

TRIPLET_DEFAULTS = {'CCA': True, 'GCA': True, 'UCA': True,
				 	'CAA': True, 'CUA': False, 'ACA': False,
					'CCU': False, 'CGA': False, 'CCC': False}
REGIONS_DEFAULTS = {'5UTR': False, 'CDS': False, '3UTR': True}
BLAST_DEFAULTS = {
	'use_blast': True,
	'db_path': str(Path.home() / 'blastdb'),
	'db_name': 'refseq_select_rna',
	'only_overlapping': False,
	'taxids': '9606',
	'taxids_placeholder': '9606 (homo sapiens)',
	'e_value': '0.05',
	'word_size': '7',
	'perc_identity': '',
	'qcov_hsp_perc': ''
}


# TODO: rename to 'validate...' etc etc or something
def on_run_analysis(controller: Controller, sequence_data: dict[str, str],
					options: dict[str, dict[str, str]]) -> None:
	# TODO: rename 'options' to something more meaningfull
	job_options = compile_job_options(options)
	job_data = {'sequence_data': sequence_data,
				'options': job_options}

	if not validate_job_options(job_options):
		return

	if not job_options['blast']['use_blast']:
		controller.start_job(source='/options',
							 job_fn=generate_sensors_job,
							 job_data=job_data)
		return

	blast_db_path = job_options['blast']['db_path']
	blast_db_name = job_options['blast']['db_name']
	blastdbcheck_binary = controller.detect_blastdbcheck_installation()

	if not blastdbcheck_binary:
		logging.warning('blastdbcheck executable not found')
		show_dialog_box(f'could not validate blast database: blastdbcheck missing')
		return

	if not os.path.exists(blast_db_path) or not os.path.isdir(blast_db_path):
		show_dialog_box(f'directory {blast_db_path} not found')
		return

	if not controller.check_blast_db(blastdbcheck_binary=blastdbcheck_binary,
									 db_path=blast_db_path,
									 db_name=blast_db_name):
		show_dialog_box(f"invalid blast database '{blast_db_name}'")
		return

	controller.start_job(source='/options',
						 job_fn=generate_sensors_job,
						 job_data=job_data)


def compile_job_options(options: dict[str, dict[str, str]]) -> dict:
	# TODO: rename 'options' to something more meaningfull
	job_options = {
		'triplets': tuple(triplet for triplet, checked in options['triplets'].items() if checked),
		'regions': tuple(region for region, checked in options['regions'].items() if checked),
		'blast': options['blast']
	}
	return job_options


def validate_job_options(job_options: dict) -> bool:
	options_to_validate: list[tuple[str, Callable, str]] = [
		(job_options['triplets'], any,
		 'at least one triplet required'),
		(job_options['regions'], any,
		 'at least one region required')
	]

	if job_options['blast']['use_blast']:
		options_to_validate.extend([
			(job_options['blast']['evalue'],
			 lambda value: float(value) > 0,
			 'E-value must be a positive number'),
			(job_options['blast']['word_size'],
			 lambda value: int(value) >= 4,
			 'word size must be an integer number >= 4'),
			(job_options['blast']['perc_identity'],
			 lambda value: True if not value else 0 <= float(value) <= 100,
			 'Min. percent identity must be a number between 0 and 100'),
			(job_options['blast']['qcov_hsp_perc'],
			 lambda value: True if not value else 0 <= float(value) <= 100,
			 'Min. percent query coverage must be a number between 0 and 100'),
		])

	for option_value, validate_fn, err_msg in options_to_validate:
		try:
			success = validate_fn(option_value)
		except ValueError:
			success = False
		if not success:
			show_dialog_box(err_msg)
			return False
	return True


@page_builder(allow_from='/metainf')
def build(controller: Controller, data: dict[str, Any] | None) -> None:
	if not data.get('sequence_data'):
		show_page_not_available(controller)
		return

	back_button(controller, route='/metainf', source='/options', data=data)

	# TODO: rename 'options' to something more meaningfull
	options = {'triplets': {},
			   'regions': {},
			   'blast': {}}

	blast_binary = controller.detect_blastn_installation()
	if not blast_binary:
		options['blast']['use_blast'] = False
	else:
		options['blast']['blast_binary'] = blast_binary

	with ui.row().classes('w-full justify-center'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-6 w-5/6 h-auto w-full'):
			with ui.row().classes('w-full place-content-center justify-around'):

				with ui.column().classes('place-content-center'):
					ui.label('Target triplets')\
						.classes('text-center text-lg font-mono font-semibold')
					for triplet in TRIPLET_DEFAULTS.keys():
						checkbox = ui.checkbox(triplet).classes('text-base font-mono')
						checkbox.value = TRIPLET_DEFAULTS[triplet]
						checkbox.bind_value(options['triplets'], triplet)

				with ui.column().classes('place-content-center'):
					ui.label('Target regions')\
						.classes('text-center text-lg font-mono font-semibold')
					for region in REGIONS_DEFAULTS.keys():
						checkbox = ui.checkbox(region).classes('text-base font-mono')
						checkbox.value = REGIONS_DEFAULTS[region]
						checkbox.bind_value(options['regions'], region)

				blast_options_column = ui.column()
				with blast_options_column.classes('place-content-center'):
					ui.label('BLAST').classes('text-center text-lg font-mono font-semibold')

					use_blast = ui.checkbox('run blast').classes('text-base font-mono')
					use_blast.value = BLAST_DEFAULTS['use_blast']
					use_blast.bind_value(options['blast'], 'use_blast')

					only_include_overlapping = ui.checkbox('only include hits overlapping with'
														   ' target triplet').classes('text-base font-mono')
					only_include_overlapping.value = BLAST_DEFAULTS['only_overlapping']
					only_include_overlapping.bind_value(options['blast'], 'only_overlapping')

					with ui.grid(columns=2):
						ui.label('word size') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						word_size_input = ui.input(value=BLAST_DEFAULTS['word_size'],
												   placeholder=BLAST_DEFAULTS['word_size'])\
							.classes('text-center text-base font-mono font-semibold')
						word_size_input.bind_value(options['blast'], 'word_size')

						ui.label('E-value threshold') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						e_value_input = ui.input(value=BLAST_DEFAULTS['e_value'],
												 placeholder=BLAST_DEFAULTS['e_value'])\
							.classes('text-center text-base font-mono font-semibold')
						e_value_input.bind_value(options['blast'], 'evalue')

						ui.label('Min. percent identity') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						perc_identity_input = ui.input(value=BLAST_DEFAULTS['perc_identity'],
												placeholder=BLAST_DEFAULTS['perc_identity'])\
							.classes('text-center text-base font-mono font-semibold')
						perc_identity_input.bind_value(options['blast'], 'perc_identity')

						ui.label('Min. percent query coverage') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						qcov_hsp_perc_input = ui.input(value=BLAST_DEFAULTS['qcov_hsp_perc'],
												placeholder=BLAST_DEFAULTS['qcov_hsp_perc'])\
							.classes('text-center text-base font-mono font-semibold')
						qcov_hsp_perc_input.bind_value(options['blast'], 'qcov_hsp_perc')

						ui.label('Filter restults by taxonomic ID:') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						taxids_input = ui.input(value=BLAST_DEFAULTS['taxids'],
												placeholder=BLAST_DEFAULTS['taxids_placeholder'])\
							.classes('text-center text-base font-mono font-semibold')
						taxids_input.bind_value(options['blast'], 'taxids')

						ui.label('Path to BLAST database directory') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						db_path_input = ui.input(value=BLAST_DEFAULTS['db_path'])\
							.classes('text-center text-base font-mono font-semibold')
						db_path_input.bind_value(options['blast'], 'db_path')

						ui.label('BLAST database name') \
							.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
						db_name_input = ui.input(value=BLAST_DEFAULTS['db_name'])\
							.classes('text-center text-base font-mono font-semibold')
						db_name_input.bind_value(options['blast'], 'db_name')

			if not blast_binary:
				blast_options_column.clear()
				with blast_options_column:
					ui.label('BLAST').classes('text-center text-lg font-mono font-semibold')
					ui.label('no blast installation detected') \
						.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')\
						.style(f'color: {Colors.RED}')

			ui.button('run analysis',
					  on_click=lambda: on_run_analysis(controller, data['sequence_data'], options))\
				.classes('w-2/3 self-center h-10 font-mono mt-6')

from nicegui import ui
from nicegui.events import UploadEventArguments
from typing import Any
from interface import Controller
from .style import Colors
from .pagebuilder import page_builder
from .ui_elements import show_dialog_box, back_button
import seqinput


MAX_INPUT_LEN: int = 500_000  # 500 kb (for comparison, the human titin TTN-018 transcript is ~109 kb long)


def show_file_size_error() -> None:
	show_dialog_box((
		'This file is too large. '
		f'Maximum file size: {int(MAX_INPUT_LEN / 1000)} KB'
	))


def show_file_decode_error() -> None:
	show_dialog_box((
		'This file could not be read. '
		'Please upload a file in FASTA or Genbank format.'
	))


def show_invalid_format_error() -> None:
	show_dialog_box((
		'This is not a single-record FASTA/GenBank sequence. '
		'Please upload a file in FASTA/Genbank format or '
		'manually enter a sequence.'
	))


def show_empty_sequence_error() -> None:
	show_dialog_box('This sequence record is empty.')


def process_file_upload(upload: UploadEventArguments, controller: Controller) -> None:
	try:
		sequence: str = upload.content.read().decode('UTF-8')
		process_sequence_input(sequence, controller)
	except UnicodeDecodeError:
		show_file_decode_error()


def process_sequence_input(sequence: str, controller: Controller) -> None:
	try:
		sequence_data: dict[str, str] = seqinput.read_sequence(sequence)
		if not sequence_data['sequence']:
			show_empty_sequence_error()
		else:
			controller.open_page(route='/metainf', source='/input',
								 data={'sequence_data': sequence_data})
	except ValueError:
		# ValueError is raised if Bio.SeqIO.read can't read the sequence as (single-entry) FASTA/GenBank record
		# or if the sequence contains letters other than A,C,G,U after transcribing to RNA
		show_invalid_format_error()


@page_builder()
def build(controller: Controller, data: dict[str, Any] | None) -> None:
	back_button(controller, '/')

	with ui.row().classes('w-full justify-center mt-24'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/5 h-64'):
			with ui.column().classes('w-full'):
				with ui.row().classes('w-full justify-center'):
					ui.label('Upload a file').classes('text-center text-lg font-mono')
					ui.image('/assets/img-file.png').classes('w-7 h-7')
				ui.upload(on_upload=lambda upload: process_file_upload(upload, controller),
						  on_rejected=show_file_size_error,
						  max_file_size=MAX_INPUT_LEN,
						  max_files=1).classes('w-full mt-10')

		ui.label('or')\
			.classes('text-2xl font-semibold self-center pl-6 pr-6 font-mono')\
			.style(f'color: {Colors.ACCENT}')

		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/5 h-64'):
			with ui.column().classes('w-full'):
				with ui.row().classes('w-full justify-center'):
					ui.label('Enter a sequence').classes('text-center text-lg font-mono')
					ui.image('/assets/img-clipboard.png').classes('w-7 h-7')

				sequence_input = ui.textarea(label='FASTA/GenBank record',
											 placeholder='>NM_005228.5 EGFR...',
											 validation={
												 'Input too long': lambda value: len(value) < MAX_INPUT_LEN})\
					.classes('w-full')

				ui.button('continue',
						  on_click=lambda: process_sequence_input(sequence_input.value, controller))\
					.classes('w-2/3 self-center h-10 font-mono')
	ui.query('textarea').style('resize: none').classes('h-24')

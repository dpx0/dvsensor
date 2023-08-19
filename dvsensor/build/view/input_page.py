import sys
sys.path.append('..')

from nicegui import ui, events
from .base_elements import header, footer, back_button
from .style import Colors, set_colors
import utils


MAX_INPUT_LEN: int = 200_000  # 200 kb (for comparison, the human titin TTN-018 transcript is ~109 kb long)


def handle_file_upload(upload: events.UploadEventArguments, view) -> None:
	try:
		content: str = upload.content.read().decode('UTF-8')
	except UnicodeDecodeError:
		show_file_decode_error()
		return
	seq = utils.read_fasta_string(content)
	if seq is not None:
		view.controller.process_fasta_sequence(seq)
	else:
		show_invalid_fasta_error()


def handle_manual_input(input: str, view) -> None:
	seq = utils.read_fasta_string(input)
	if seq is not None:
		view.controller.process_fasta_sequence(seq)
	else:
		show_invalid_fasta_error()


def show_file_size_error() -> None:
	with ui.dialog() as dialog, ui.card():
		ui.label('This file is too large').classes('text-lg self-center font-mono')
		ui.label(f'Maximum file size: {int(MAX_INPUT_LEN/1000)} KB').classes('text-lg self-center font-mono')
		ui.button('OK', on_click=dialog.close).classes('w-full self-center text-sm font-mono')
	dialog.open()


def show_file_decode_error() -> None:
	with ui.dialog() as dialog, ui.card():
		ui.label('This file could not be read').classes('text-lg self-center font-mono')
		ui.label(f'Please upload a file in FASTA format').classes('text-lg self-center font-mono')
		ui.button('OK', on_click=dialog.close).classes('w-full self-center text-sm font-mono')
	dialog.open()


def show_invalid_fasta_error() -> None:
	with ui.dialog() as dialog, ui.card():
		ui.label('This is not a FASTA sequence').classes('text-lg self-center font-mono')
		ui.label('Please upload a file in FASTA format or').classes('text-lg self-center font-mono')
		ui.label('manually enter a equence').classes('text-lg self-center font-mono')
		ui.button('OK', on_click=dialog.close).classes('w-full self-center text-sm font-mono')
	dialog.open()


def build(view, **kwargs) -> None:
	set_colors()
	header()
	back_button('start', view)

	with ui.row().classes('w-full justify-center mt-24'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/5 h-64'):
			with ui.column().classes('w-full'):
				with ui.row().classes('w-full justify-center'):
					ui.label('Upload a FASTA file').classes('text-center text-lg font-mono')
					ui.image('/assets/img-file.png').classes('w-7 h-7')
				ui.upload(on_upload=lambda upload: handle_file_upload(upload, view),
						  on_rejected=show_file_size_error,
						  max_file_size=MAX_INPUT_LEN,
						  max_files=1).classes('w-full mt-10')

		ui.label("or"
				 ).classes('text-2xl font-semibold self-center pl-6 pr-6 font-mono'
				 ).style(f'color: {Colors.ACCENT}')

		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/5 h-64'):
			with ui.column().classes('w-full'):
				with ui.row().classes('w-full justify-center'):
					ui.label('Enter a FASTA sequence').classes('text-center text-lg font-mono')
					ui.image('/assets/img-clipboard.png').classes('w-7 h-7')
				sequence_input = ui.textarea(label='FASTA Sequence', placeholder='>NM_005228.5 EGFR...',
										  validation={'Input too long': lambda value: len(value) < MAX_INPUT_LEN}
										  ).classes('w-full')
				ui.button('continue', on_click=lambda: handle_manual_input(sequence_input.value, view)
						  ).classes('w-2/3 self-center h-10 font-mono')
	footer()




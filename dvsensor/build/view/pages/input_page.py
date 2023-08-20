from nicegui import ui
from ..base_elements import header, footer, back_button, show_dialog_box
from ..style import Colors, set_colors


MAX_INPUT_LEN: int = 200_000  # 200 kb (for comparison, the human titin TTN-018 transcript is ~109 kb long)


def handle_seq_file_upload(upload, view) -> None:
	try:
		user_input: str = upload.content.read().decode('UTF-8')
		view.controller.handle_fasta_seq_input(user_input)
	except UnicodeDecodeError:
		show_file_decode_error()
	except ValueError:  # raised when Bio.SeqIO.read(...) can't parse a string as a FASTA record
		show_invalid_fasta_error()


def handle_manual_seq_input(user_input: str, view) -> None:
	try:
		view.controller.handle_fasta_seq_input(user_input)
	except ValueError:  # raised when Bio.SeqIO.read(...) can't parse a string as a FASTA record
		show_invalid_fasta_error()


def show_file_size_error() -> None:
	show_dialog_box([
		'This file is too large',
		f'Maximum file size: {int(MAX_INPUT_LEN / 1000)} KB'
	])


def show_file_decode_error() -> None:
	show_dialog_box([
		'This file could not be read',
		'Please upload a file in FASTA format'
	])


def show_invalid_fasta_error() -> None:
	show_dialog_box([
		'This is not a single-record', 'FASTA sequence',
		'Please upload a file in FASTA format or',
		'manually enter a sequence'
	])


def build(view, **kwargs) -> None:
	view.controller.clear_job()
	set_colors()
	header()
	back_button('start', view)

	with ui.row().classes('w-full justify-center mt-24'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/5 h-64'):
			with ui.column().classes('w-full'):
				with ui.row().classes('w-full justify-center'):
					ui.label('Upload a FASTA file').classes('text-center text-lg font-mono')
					ui.image('/assets/img-file.png').classes('w-7 h-7')
				ui.upload(on_upload=lambda upload: handle_seq_file_upload(upload, view),
						  on_rejected=show_file_size_error,
						  max_file_size=MAX_INPUT_LEN,
						  max_files=1).classes('w-full mt-10')

		ui.label('or')\
			.classes('text-2xl font-semibold self-center pl-6 pr-6 font-mono')\
			.style(f'color: {Colors.ACCENT}')

		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/5 h-64'):
			with ui.column().classes('w-full'):
				with ui.row().classes('w-full justify-center'):
					ui.label('Enter a FASTA sequence').classes('text-center text-lg font-mono')
					ui.image('/assets/img-clipboard.png').classes('w-7 h-7')
				sequence_input = ui.textarea(label='FASTA Sequence', placeholder='>NM_005228.5 EGFR...',
											 validation={'Input too long': lambda value: len(value) < MAX_INPUT_LEN})\
					.classes('w-full')
				ui.button('continue',
						  on_click=lambda: handle_manual_seq_input(sequence_input.value, view))\
					.classes('w-2/3 self-center h-10 font-mono')

	footer()

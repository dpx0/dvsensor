from nicegui import ui, events
from .base_elements import header, footer, back_button
from .style import Colors, set_colors


MAX_INPUT_LEN: int = 200_000  # 200 kb (for comparison, the human titin TTN-018 transcript is ~109 kb long)


def handle_file_upload(upload: events.UploadEventArguments, view) -> None:
	view.controller.process_input_sequence(upload.content.read().decode('UTF-8'))


def handle_file_too_large() -> None:
	with ui.dialog() as dialog, ui.card():
		ui.label('This file is too large!').classes('text-lg self-center font-mono')
		ui.label(f'Maximum file size: {int(MAX_INPUT_LEN/1000)} KB').classes('text-lg self-center font-mono')
		ui.button('OK', on_click=dialog.close).classes('w-full self-center text-sm font-mono')
	dialog.open()


def handle_manual_input(sequence: str, view) -> None:
	view.controller.process_input_sequence(sequence)


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
						  on_rejected=handle_file_too_large,
						  max_file_size=MAX_INPUT_LEN).classes('w-full mt-10')

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




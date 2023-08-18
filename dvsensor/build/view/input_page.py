from nicegui import ui, events


def handle_file_upload(upload: events.UploadEventArguments, view) -> None:
	view.controller.process_input_sequence(upload.content.read().decode('UTF-8'))


def build(view, **kwargs) -> None:
	ui.label('Welcome to the Input Page')
	ui.label(f'View reference {id(view)}')
	ui.label('Upload a FASTA file:')
	ui.upload(on_upload=lambda upload: handle_file_upload(upload, view)).classes('max-w-full')
	ui.label('or manually enter a FASTA Sequence:')
	manual_input = ui.input(label='FASTA Sequence', placeholder='FASTA Sequence',
         validation={'Input too long': lambda value: len(value) < 20})
	ui.button('Next',
			  on_click=lambda: view.controller.process_input_sequence(manual_input.value)).classes('w-32')



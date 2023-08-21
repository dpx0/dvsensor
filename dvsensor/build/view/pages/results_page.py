from nicegui import ui
from ..base_elements import header, footer
from ..style import Colors, set_colors
from functools import partial


def add_rows(grid, row_data: list[dict]) -> None:
	grid.call_api_method('applyTransaction', {'add': row_data})


def update_progress(progress_bar, progress_label, step: float) -> float:
	progress_bar.value += step
	progress_label.text = f'{progress_bar.value * 100: .0f} %'
	return progress_bar.value


def set_status_finished(spinner, status_text):
	spinner.visible = False
	status_text.text = 'Analysis finished'
	status_text.style(replace=f'color: {Colors.GREEN}')


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if not view.controller.page_allowed(kwargs):
		view.show_error('page not available')
		return

	with ui.column().classes('w-full').style('height: 80vh'):

		# ----- top row
		with ui.row().classes('w-full flex justify-between mt-2 mb-10 space-x-10'):

			# ----- target ideogram
			with ui.column().classes('grow'):
				ui.label(f'Target: {view.model.record_name} ({view.model.record_id})')\
					.classes('text-center text-lg font-mono font-semibold')\
					.style(f'color: {Colors.FOREGROUND}')

			# ----- target information
			with ui.column().classes('1/4'):
				ui.label(f'Total transcript length: {view.model.sequence_length} bp')\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				ui.label(f"5'-UTR length: {view.model.sequence_length} bp")\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				ui.label(f'CDS length: {view.model.sequence_length} bp')\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				ui.label(f"3'-UTR length: {view.model.sequence_length} bp")\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')

			# ----- status + control elements
			with ui.column().classes('w-1/3 flex justify-between place-content-center'):
				with ui.row():
					spinner = ui.spinner(color=Colors.ACCENT, size='xl').classes()
					status_text = ui.label(f'Analysis running...')\
						.classes('self-center text-lg font-mono font-semibold')\
						.style(f'color: {Colors.ACCENT}')

		# ----- progress bar + label
		with ui.grid(columns=2).classes('w-full flex space-x-4 pl-6'):
			progress_label = ui.label('0 %')\
				.classes('w-16 text-lg font-mono font-semibold')\
				.style(f'color: {Colors.FOREGROUND}')
			progress_bar = ui.linear_progress(show_value=False)\
				.props('track-color=grey-1 color=green-9')\
				.classes('w-0 h-4 grow self-center')

		# ----- sensor table
		sensors_table = ui.aggrid({
			'defaultColDef': {'flex': 1},
			'columnDefs': [
				{'headerName': 'Position', 'field': 'position'},
				{'headerName': 'Triplet', 'field': 'triplet'},
				{'headerName': 'Region', 'field': 'region'},
				{'headerName': 'Range', 'field': 'range'},
				{'headerName': '%GC', 'field': 'percent_gc'},
				{'headerName': '# In-frame stop codons', 'field': 'n_stop_codons'},
				{'headerName': 'Off-targets', 'field': 'off_targets'},
			],
			'rowData': [],
			'rowSelection': 'multiple',
		}, theme='alpine-dark').classes('w-full grow')

		task_control_functions = view.controller.task_controller.start_analysis_task(
			ui_control_functions={
				'add_rows': partial(add_rows, sensors_table),
				'update_progress': partial(update_progress, progress_bar, progress_label),
				'set_status_finished': lambda: set_status_finished(spinner, status_text)
			})

		# ui.button('CANCEL', on_click=task_control_functions['cancel_task'])

	footer()

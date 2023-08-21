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


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if not view.controller.page_allowed(kwargs):
		view.show_error('page not available')
		return

	with ui.column().classes('w-full flex').style('height: 80vh'):

		with ui.row().classes('w-full justify-between mt-2 mb-10'):
			with ui.column().classes('ml-6'):
				ui.label(f'Target: {view.model.record_name} ({view.model.record_id})')\
					.classes('text-center text-lg font-mono font-semibold')\
					.style(f'color: {Colors.FOREGROUND}')

		with ui.grid(columns=2).classes('w-full flex pl-6 space-x-4'):
			progress_label = ui.label('0 %')\
				.classes('w-16 text-lg font-mono font-semibold')\
				.style(f'color: {Colors.FOREGROUND}')
			progress_bar = ui.linear_progress(show_value=False)\
				.props('track-color=grey-1 color=green-9')\
				.classes('w-0 h-4 grow self-center')

		sensors_table = ui.aggrid({
			'defaultColDef': {'flex': 1},
			'columnDefs': [
				{'headerName': 'Position', 'field': 'position'},
				{'headerName': 'Triplet', 'field': 'triplet'},
				{'headerName': 'Region', 'field': 'region'},
				{'headerName': 'Range', 'field': 'range'},
				{'headerName': '%GC', 'field': 'percent_gc'},
				{'headerName': '# Stop Codon Edits', 'field': 'n_stop_edits'},
				{'headerName': 'off-targets', 'field': 'off_targets'},
			],
			'rowData': [],
			'rowSelection': 'multiple',
		}, theme='alpine-dark').classes('w-full grow')

		job_control_functions = view.controller.start_analysis_job(ui_control_functions={
			'add_rows': partial(add_rows, sensors_table),
			'update_progress': partial(update_progress, progress_bar, progress_label)
		})

		ui.button('CANCEL', on_click=job_control_functions['cancel_job'])

	footer()

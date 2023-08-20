from nicegui import ui
from ..base_elements import header, footer
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()
	#if not view.model.has_record_data:
		#view.show_error('page not available')
		#return

	with ui.column().classes('w-full'):

		with ui.row().classes('w-full justify-between mt-2 mb-10'):
			with ui.column().classes('ml-6'):
				ui.label(f'Target: {view.model.record_name} ({view.model.record_id})')\
					.classes('text-center text-lg font-mono font-semibold')\
					.style(f'color: {Colors.FOREGROUND}')

		with ui.grid(columns=2).classes('w-full flex pl-6 space-x-6'):
			progress_label = ui.label('100%')\
				.classes('w-5 text-lg font-mono font-semibold')\
				.style(f'color: {Colors.FOREGROUND}')
			progress_bar = ui.linear_progress(show_value=False)\
				.props('track-color=grey-1 color=green-9')\
				.classes('w-0 h-4 grow self-center')

		grid = ui.aggrid({
			'defaultColDef': {'flex': 1},
			'columnDefs': [
				{'headerName': 'Position', 'field': 'position'},
				{'headerName': 'Triplet', 'field': 'triplet'},
				{'headerName': 'Region', 'field': 'region'},
				{'headerName': 'Range', 'field': 'range'},
				{'headerName': '%GC', 'field': 'percent_gc'},
				{'headerName': '# Stop Codon Edits', 'field': 'n_stop_edits'},
				{'headerName': 'off-targets', 'field': 'off_targets'},
				{'headerName': '', 'field': 'buttons'}
			],
			'rowData': [
				{'name': 'Alice', 'age': 18, 'parent': 'David'},
				{'name': 'Bob', 'age': 21, 'parent': 'Eve'},
				{'name': 'Carol', 'age': 42, 'parent': 'Frank'},
			],
			'rowSelection': 'multiple',
		}, theme='alpine-dark').classes('w-full')

	footer()

from nicegui import ui
from functools import partial
from typing import Any
from interface import Controller
from signals import Signal
from .pagebuilder import page_builder
from .ui_elements import show_page_not_available, show_dialog_box
from .style import Colors


def add_rows(output_table: ui.aggrid, rows: list[dict]) -> None:
	output_table.call_api_method('applyTransaction', {'add': rows})


def update_progress(progress_bar: ui.linear_progress, progress_label: ui.label, step: float) -> None:
	progress_bar.value += step
	if round(progress_bar.value * 100) >= 100:
		progress_bar.value = 1.0
	progress_label.text = f'{progress_bar.value * 100.0: .0f} %'


def set_status_finished(spinner: ui.spinner, cancel_button: ui.button,
						export_button: ui.button, status_label: ui.label) -> None:
	print("set status finished")
	spinner.visible = False
	cancel_button.visible = False
	export_button.visible = True
	status_label.text = 'Analysis finished'
	status_label.style(replace=f'color: {Colors.GREEN}')


def set_status_cancelled(spinner: ui.spinner, cancel_button: ui.button,
						export_button: ui.button, status_label: ui.label) -> None:
	spinner.visible = False
	cancel_button.visible = False
	export_button.visible = True
	status_label.text = 'Analysis cancelled'
	status_label.style(replace=f'color: {Colors.RED}')


def show_confirm_cancel_job_dialog(controller: Controller, spinner: ui.spinner, cancel_button: ui.button,
								   export_button: ui.button, status_label: ui.label) -> None:
	buttons = show_dialog_box('Are you sure you want to cancel the analysis?',
							  ['yes', 'no'])
	buttons['yes'].on('click', lambda: cancel_job(controller, spinner, cancel_button,
												  export_button, status_label))


def cancel_job(controller: Controller, spinner: ui.spinner, cancel_button: ui.button,
			   export_button: ui.button, status_label: ui.label) -> None:
	set_status_cancelled(spinner, cancel_button, export_button, status_label)
	controller.terminate_job()


def show_confirm_home_dialog(controller: Controller, spinner: ui.spinner, cancel_button: ui.button,
							 export_button: ui.button, status_label: ui.label) -> None:
	if not controller.is_job_running():
		controller.open_page('/')
	else:
		buttons = show_dialog_box((
			'Are you sure you want to leave the page? '
			'This will cancel the running analysis.')
			, ['yes', 'no'])
		buttons['yes'].on('click', lambda: (cancel_job(controller, spinner, cancel_button,
													  export_button, status_label),
										   controller.open_page('/')))


async def update_ideogram_sensor_window(box_5UTR: ui.column,
										box_3UTR: ui.column,
										sensor_range_str: str,
										transcript_len: int) -> None:
	sensor_start, sensor_end = sensor_range_str.split('-')
	sensor_start = int(sensor_start)
	sensor_end = int(sensor_end)
	await ui.run_javascript(f'''
	var box_5UTR = getElement({box_5UTR.id});
	var box_3UTR = getElement({box_3UTR.id});
	
	var offsets_box_5UTR = box_5UTR.getBoundingClientRect();
	var offsets_box_3UTR = box_3UTR.getBoundingClientRect();
	
	var ideog_top = offsets_box_5UTR.top + document.documentElement.scrollTop;
	var ideog_bottom = offsets_box_5UTR.bottom + document.documentElement.scrollTop;
	var ideog_left = offsets_box_5UTR.left + document.documentElement.scrollLeft;
	var ideog_right = offsets_box_3UTR.right + document.documentElement.scrollLeft;
	var ideog_len = ideog_right - ideog_left;

	if (document.getElementById("sensor_window") == undefined) {{	
		var sensorWindow = document.createElement("div");
		sensorWindow.id = "sensor_window"
		sensorWindow.className = "non-ng"

		sensorWindow.style = "background-color: rgba(239,68,68,.4)";
		sensorWindow.style.border = "medium solid #ef4444"
		sensorWindow.style.position = "absolute";
		sensorWindow.style.top = ideog_top + 'px';
		sensorWindow.style.height = ideog_bottom - ideog_top + 'px';

		document.body.appendChild(sensorWindow);
	}} else {{
		var sensorWindow = document.getElementById("sensor_window");
	}}

	sensorWindow.style.left = ideog_left + (ideog_len * ({sensor_start / transcript_len})) + 'px';
	sensorWindow.style.width =  (ideog_len * ({(sensor_end - sensor_start) / transcript_len})) + 'px';
	''')


@page_builder(allow_from='/options')
def build(controller: Controller, data: dict[str, Any] | None) -> None:
	if not data.get('sequence_data'):
		show_page_not_available(controller)
		return

	ideogram_data_ready_signal: Signal = controller.signal_handler.request_signal('ideogram_data_ready')
	ideogram_data: dict[str, int | float] = ideogram_data_ready_signal.wait()
	len_5UTR = ideogram_data['len_5UTR']
	len_CDS = ideogram_data['len_CDS']
	len_3UTR = ideogram_data['len_3UTR']
	percent_5UTR = ideogram_data['percent_5UTR']
	percent_CDS = ideogram_data['percent_CDS']
	percent_3UTR = ideogram_data['percent_3UTR']

	controller.signal_handler.register_signal('ui_ready')

	with ui.column().classes('w-full').style('height: 80vh'):

		# ----- top row
		with ui.row().classes('w-full flex justify-between mt-2 mb-10 space-x-10'):

			# ----- target name + ideogram
			with ui.column().classes('grow'):

				# ----- target name
				ui.label(f'Target: {data["sequence_data"]["name"]}' +
						 f'({data["sequence_data"]["accession"]})')\
					.classes('text-center text-lg font-mono font-semibold')\
					.style(f'color: {Colors.FOREGROUND}')

				# ----- target ideogram
				with ui.row().classes('w-full gap-0 pt-8'):
					box_5UTR = ui.column().style(f'width: {percent_5UTR}%')
					with box_5UTR:
						with ui.element('div').classes('w-full bg-stone-500'):
							ui.label(f"""{"5'-UTR" if percent_5UTR >= 10 else '*'}""") \
								.classes('text-center font-semibold')

					box_CDS = ui.column().style(f'width: {percent_CDS}%')
					with box_CDS:
						with ui.element('div').classes('w-full bg-yellow-500'):
							ui.label(f"{'CDS' if percent_CDS >= 10 else '*'}") \
								.classes('text-center font-semibold')

					box_3UTR = ui.column().style(f'width: {percent_3UTR}%')
					with box_3UTR:
						with ui.element('div').classes('w-full bg-indigo-500'):
							ui.label(f"""{"3'-UTR" if percent_3UTR >= 10 else '*'}""") \
								.classes('text-center font-semibold')

				with ui.row().classes('w-full justify-center gap-x-28'):
					ui.label(f"5'-UTR: {len_5UTR} bp") \
						.classes('text-center text-sm font-mono') \
						.style(f'color: {Colors.FOREGROUND}')
					ui.label(f"CDS: {len_CDS} bp") \
						.classes('text-center text-sm font-mono') \
						.style(f'color: {Colors.FOREGROUND}')
					ui.label(f"3'-UTR: {len_3UTR} bp") \
						.classes('text-center text-sm font-mono') \
						.style(f'color: {Colors.FOREGROUND}')

			# ----- additional information
			with ui.column().classes('1/4'):
				ui.label(f'Total transcript length: {len_3UTR + len_CDS + len_5UTR} bp')\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				ui.label(f'Target triplets analyzed:')\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')

			# ----- status + control elements
			status_box = ui.column().classes('w-1/3 flex justify-between place-content-center')
			with status_box:
				with ui.row():
					spinner = ui.spinner(color=Colors.ACCENT, size='xl').classes()
					status_label = ui.label(f'Analysis running...')\
						.classes('self-center text-lg font-mono font-semibold')\
						.style(f'color: {Colors.ACCENT}')

		ui.separator().props('dark')

		# ----- progress bar + label
		with ui.grid(columns=2).classes('w-full flex space-x-4 pl-6'):
			progress_label = ui.label('0 %')\
				.classes('w-16 text-lg font-mono font-semibold')\
				.style(f'color: {Colors.FOREGROUND}')
			progress_bar = ui.linear_progress(show_value=False)\
				.props('track-color=grey-1 color=green-9 stripe')\
				.classes('w-0 h-4 grow self-center')

		# ----- sensor table
		output_table = ui.aggrid({
			'defaultColDef': {'flex': 1},
			'columnDefs': [
				{'headerName': 'Position', 'field': 'position'},
				{'headerName': 'Triplet', 'field': 'triplet'},
				{'headerName': 'Region', 'field': 'region'},
				{'headerName': 'Range', 'field': 'range'},
				{'headerName': '%GC', 'field': 'percent_gc'},
				{'headerName': 'In-frame start/stop codons', 'field': 'n_edits'},
				{'headerName': 'Potential off-targets', 'field': 'pot_off_targets'},
				{'headerName': 'Sensor (5->3)', 'field': 'sensor'},
				{'headerName': 'Trigger (5->3)', 'field': 'trigger'}
			],
			'rowData': [],
			'rowSelection': 'multiple',
		}, theme='alpine-dark')\
			.classes('w-full grow')\
			.on('rowSelected',
				 lambda event: update_ideogram_sensor_window(
					 box_5UTR, box_3UTR,
					 event.args['data']['range'],
					 len_3UTR + len_CDS + len_5UTR)
				if event.args['selected'] else None)

		with status_box:
			export_button = ui.button('export table as CSV',
									  on_click=lambda: output_table.call_api_method(
										  'exportDataAsCsv', {'suppressQuotes': True})) \
				.props('color=indigo-10')\
				.classes('self-center')
			export_button.visible = False

			cancel_button = ui.button('cancel',
									  on_click=lambda: show_confirm_cancel_job_dialog(
										  controller, spinner, cancel_button, export_button, status_label))\
				.props('color=red-10')\
				.classes('self-center')

	# ----- home button
	with ui.button('', color=f'{Colors.ACCENT}',
				   on_click=lambda: show_confirm_home_dialog(
					   controller, spinner, cancel_button, export_button, status_label))\
		.props(f'no-caps flat round size=md')\
		.classes('mt-4'):
		ui.image('/assets/img-home.png').classes('m-1')

	ui_functions = {
		'update_progress': partial(update_progress, progress_bar, progress_label),
		'set_status_finished': partial(set_status_finished, spinner, cancel_button, export_button, status_label),
		'add_rows': partial(add_rows, output_table)
	}
	controller.signal_handler.send_signal('ui_ready', args=ui_functions)

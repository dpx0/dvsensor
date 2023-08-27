from nicegui import ui
from ..base_elements import header, footer, home_button
from ..style import Colors, set_colors


def add_rows(ui_connection, row_data: list[dict]) -> None:
	ui_connection.get_element('output_table').call_api_method('applyTransaction', {'add': row_data})


def update_progress(ui_connection, step: float) -> float:
	progress_bar = ui_connection.get_element('progress_bar')
	progress_bar.value += step
	ui_connection.get_element('progress_label').text = f'{progress_bar.value * 100.0: .0f} %'


def set_status_finished(ui_connection) -> None:
	ui_connection.get_element('spinner').visible = False
	ui_connection.get_element('cancel_button').visible = False
	ui_connection.get_element('export_button').visible = True
	ui_connection.get_element('status_text').text = 'Analysis finished'
	ui_connection.get_element('status_text').style(replace=f'color: {Colors.GREEN}')


def set_status_cancelled(ui_connection) -> None:
	ui_connection.get_element('spinner').visible = False
	ui_connection.get_element('cancel_button').visible = False
	ui_connection.get_element('export_button').visible =  True
	ui_connection.get_element('status_text').text = 'Analysis cancelled'
	ui_connection.get_element('status_text').style(replace=f'color: {Colors.RED}')


def cancel_task(ui_connection):
	ui_connection.call('cancel_task')
	set_status_cancelled(ui_connection)


async def update_ideogram_sensor_window(box_5UTR, box_3UTR, sensor_range: str,
										transcript_len: int) -> None:
	sensor_start, sensor_end = sensor_range.split('-')
	sensor_start = int(sensor_start)
	sensor_end = int(sensor_end)
	await ui.run_javascript(f'''
	var offsets_box_5UTR = getElement({box_5UTR.id}).getBoundingClientRect();
	var offsets_box_3UTR = getElement({box_3UTR.id}).getBoundingClientRect();
	var ideog_top = offsets_box_5UTR.top;
	var ideog_bottom = offsets_box_5UTR.bottom;
	var ideog_left = offsets_box_5UTR.left;
	var ideog_right = offsets_box_3UTR.right;
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


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if not view.page_allowed(kwargs):
		view.show_error('page not available')
		return

	ui_connection = view.controller.connect_ui()
	ui_connection.add_function('add_rows', add_rows)
	ui_connection.add_function('update_progress', update_progress)
	ui_connection.add_function('set_status_finished', set_status_finished)

	with ui.column().classes('w-full').style('height: 80vh'):

		# ----- top row
		with ui.row().classes('w-full flex justify-between mt-2 mb-10 space-x-10'):

			# ----- target name + ideogram
			with ui.column().classes('grow'):

				# ----- target name
				ui.label(f'Target: {view.controller.query_model("record_name")}' +
						 f'({view.controller.query_model("record_id")})')\
					.classes('text-center text-lg font-mono font-semibold')\
					.style(f'color: {Colors.FOREGROUND}')

				len_5UTR = 261
				len_CDS = 3633
				len_3UTR = 6011
				len_tot = len_5UTR + len_CDS + len_3UTR

				percent_5UTR = len_5UTR / len_tot * 100
				percent_CDS = len_CDS / len_tot * 100
				percent_3UTR = 100.0 - (percent_5UTR + percent_CDS)

				# ----- target ideogram
				with ui.row().classes('w-full gap-0 pt-8'):

					box_5UTR = ui.column().classes().style(f'width: {percent_5UTR}%')
					ui_connection.add_ui_element('box_5UTR', box_5UTR)  # NEW!
					with box_5UTR:
						with ui.element('div').classes('w-full bg-stone-500'):
							ui.label(f"""{"5'-UTR" if percent_5UTR >= 10 else '*'}""")\
								.classes('text-center font-semibold')

					box_CDS = ui.column().classes().style(f'width: {percent_CDS}%')
					ui_connection.add_ui_element('box_5UTR', box_CDS)  # NEW!
					with box_CDS:
						with ui.element('div').classes('w-full bg-yellow-500'):
							ui.label(f"{'CDS' if percent_CDS >= 10 else '*'}")\
								.classes('text-center font-semibold')

					box_3UTR = ui.column().classes().style(f'width: {percent_3UTR}%')
					ui_connection.add_ui_element('box_5UTR', box_3UTR)  # NEW!
					with box_3UTR:
						with ui.element('div').classes('w-full bg-indigo-500'):
							ui.label(f"""{"3'-UTR" if percent_3UTR >= 10 else '*'}""")\
								.classes('text-center font-semibold')

				with ui.row().classes('w-full justify-center gap-x-28'):
					ui.label(f"5'-UTR: {len_5UTR} bp")\
						.classes('text-center text-sm font-mono') \
						.style(f'color: {Colors.FOREGROUND}')
					ui.label(f"CDS: {len_CDS} bp")\
						.classes('text-center text-sm font-mono')\
						.style(f'color: {Colors.FOREGROUND}')
					ui.label(f"3'-UTR: {len_3UTR} bp")\
						.classes('text-center text-sm font-mono')\
						.style(f'color: {Colors.FOREGROUND}')

			# ----- additional information
			with ui.column().classes('1/4'):
				ui.label(f'Total transcript length: {view.controller.query_model("sequence_length")} bp')\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				ui.label(f"Signal Peptide: ")\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				ui.label(f'Target sites analyzed: x / y')\
					.classes('text-center text-sm font-mono')\
					.style(f'color: {Colors.FOREGROUND}')
				#ui.label(f"3'-UTR length: {view.model.sequence_length} bp")\
					#.classes('text-center text-sm font-mono')\
					#.style(f'color: {Colors.FOREGROUND}')

			# ----- status + control elements
			status_box = ui.column().classes('w-1/3 flex justify-between place-content-center')
			with status_box:
				with ui.row():
					spinner = ui.spinner(color=Colors.ACCENT, size='xl').classes()
					ui_connection.add_ui_element('spinner', spinner)
					status_text = ui.label(f'Analysis running...')\
						.classes('self-center text-lg font-mono font-semibold')\
						.style(f'color: {Colors.ACCENT}')
					ui_connection.add_ui_element('status_text', status_text)

		ui.separator().props('dark')

		# ----- progress bar + label
		with ui.grid(columns=2).classes('w-full flex space-x-4 pl-6'):
			progress_label = ui.label('0 %')\
				.classes('w-16 text-lg font-mono font-semibold')\
				.style(f'color: {Colors.FOREGROUND}')
			ui_connection.add_ui_element('progress_label', progress_label)
			progress_bar = ui.linear_progress(show_value=False)\
				.props('track-color=grey-1 color=green-9 stripe')\
				.classes('w-0 h-4 grow self-center')
			ui_connection.add_ui_element('progress_bar', progress_bar)

		# ----- sensor table
		output_table = ui.aggrid({
			'defaultColDef': {'flex': 1},
			'columnDefs': [
				{'headerName': 'Position', 'field': 'position'},
				{'headerName': 'Triplet', 'field': 'triplet'},
				{'headerName': 'Region', 'field': 'region'},
				{'headerName': 'Range', 'field': 'range'},
				{'headerName': '%GC', 'field': 'percent_gc'},
				{'headerName': 'In-frame stop codons', 'field': 'n_stop_codons'},
				{'headerName': 'Off-targets', 'field': 'off_targets'},
			],
			'rowData': [],
			'rowSelection': 'multiple',
		}, theme='alpine-dark')\
			.classes('w-full grow')\
			.on('rowSelected',
				 lambda event: update_ideogram_sensor_window(
				 	box_5UTR, box_3UTR, event.args['data']['range'], len_tot)
				if event.args['selected'] else None)
		ui_connection.add_ui_element('output_table', output_table)

		with status_box:
			export_button = ui.button('export table as CSV',
									  on_click=lambda: output_table.call_api_method(
										  'exportDataAsCsv', {'suppressQuotes': True})) \
				.props('color=indigo-10').classes('self-center')
			export_button.visible = False
			ui_connection.add_ui_element('export_button', export_button)

			cancel_button = ui.button('cancel',
									  on_click=lambda: cancel_task(ui_connection))\
				.props('color=red-10').classes('self-center')
			ui_connection.add_ui_element('cancel_button', cancel_button)

	home_button(view)
	footer()

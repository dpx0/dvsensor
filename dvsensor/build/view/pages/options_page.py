from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if not view.page_allowed(kwargs):
		view.show_error('page not available')
		return
	back_button('metainf', view, on_click=view.controller.reload_model,
				task_id=kwargs.get('task_id'))

	triplet_default_settings = {'CCA': True, 'GCA': True, 'UCA': True,
								'CAA': True, 'CUA': False, 'ACA': False,
								'CCU': False, 'CGA': False, 'CCC': False}

	regions_default_settings = {'5UTR': False, 'CDS': False, '3UTR': True}

	with ui.row().classes('w-full justify-center'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-5/6 h-auto w-full'):
			with ui.row().classes('w-full place-content-center justify-around'):

				with ui.column().classes('place-content-center'):
					ui.label('Target triplets')\
						.classes('text-center text-lg font-mono font-semibold')

					cca_checkbox = ui.checkbox('CCA').classes('text-base font-mono')
					cca_checkbox.value = triplet_default_settings['CCA']

					gca_checkbox = ui.checkbox('GCA').classes('text-base font-mono')
					gca_checkbox.value = triplet_default_settings['GCA']

					uca_checkbox = ui.checkbox('UCA').classes('text-base font-mono')
					uca_checkbox.value = triplet_default_settings['UCA']

					caa_checkbox = ui.checkbox('CAA').classes('text-base font-mono')
					caa_checkbox.value = triplet_default_settings['CAA']

					cua_checkbox = ui.checkbox('CUA').classes('text-base font-mono')
					cua_checkbox.value = triplet_default_settings['CUA']

					aca_checkbox = ui.checkbox('ACA').classes('text-base font-mono')
					aca_checkbox.value = triplet_default_settings['ACA']

					ccu_checkbox = ui.checkbox('CCU').classes('text-base font-mono')
					ccu_checkbox.value = triplet_default_settings['CCU']

					cga_checkbox = ui.checkbox('CGA').classes('text-base font-mono')
					cga_checkbox.value = triplet_default_settings['CGA']

					ccc_checkbox = ui.checkbox('CCC').classes('text-base font-mono')
					ccc_checkbox.value = triplet_default_settings['CCC']

				with ui.column().classes('place-content-center'):
					ui.label('Target regions')\
						.classes('text-center text-lg font-mono font-semibold')

					r5utr_checkbox = ui.checkbox("5'-UTR").classes('text-base font-mono')
					r5utr_checkbox.value = regions_default_settings['5UTR']

					rCDS_checkbox = ui.checkbox("CDS").classes('text-base font-mono')
					rCDS_checkbox.value = regions_default_settings['CDS']

					r3utr_checkbox = ui.checkbox("3'-UTR").classes('text-base font-mono')
					r3utr_checkbox.value = regions_default_settings['3UTR']

				with ui.column().classes('place-content-center'):
					ui.label('BLAST').classes('text-center text-lg font-mono font-semibold')

					# TODO: implement options for BLAST

				with ui.column().classes('place-content-center'):
					ui.label('Signal Peptide').classes('text-center text-lg font-mono font-semibold')

					# TODO: implement options for signal peptide prediction

			ui.button('run analysis',
					  on_click=lambda: view.controller.start_task({
						  'type': 'single',
						  'triplet_settings': {
							  'CCA': cca_checkbox.value,
							  'GCA': gca_checkbox.value,
							  'UCA': uca_checkbox.value,
							  'CAA': caa_checkbox.value,
							  'CUA': cua_checkbox.value,
							  'ACA': aca_checkbox.value,
							  'CCU': ccu_checkbox.value,
							  'CGA': cga_checkbox.value,
							  'CCC': ccc_checkbox.value
						  },
						  'regions_settings': {
							  '5UTR': r5utr_checkbox.value,
							  'CDS': rCDS_checkbox.value,
							  '3UTR': r3utr_checkbox.value
						  }
			})).classes('w-2/3 self-center h-10 font-mono mt-6')

	footer()

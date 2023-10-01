from nicegui import ui
from typing import Any
from interface import Controller
from .pagebuilder import page_builder
from .ui_elements import back_button, show_page_not_available


@page_builder(allow_from=['/input', '/options'])
def build(controller: Controller, data: dict[str, Any] | None) -> None:
	if not data.get('sequence_data'):
		show_page_not_available(controller)
		return

	back_button(controller, '/input')

	with ui.row().classes('w-full justify-center mt-24'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/3 h-auto '):
			with ui.column().classes('w-full place-content-center'):

				ui.label('Please check if the sequence record information is correct')\
					.classes('text-center text-lg font-mono')

				with ui.grid(columns=2).classes('place-self-start w-2/3'):
					ui.label('Name')\
						.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
					name_input = ui.input().classes('text-center text-base font-mono font-semibold')
					name_input.value = data['sequence_data'].get('name', '')
					name_input.bind_value(data['sequence_data'], 'name')

					ui.label('Accession Number')\
						.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
					accession_input = ui.input().classes('text-center text-base font-mono font-semibold')
					accession_input.value = data['sequence_data'].get('accession', '')
					accession_input.bind_value(data['sequence_data'], 'accession')

				ui.button('ok',
						  on_click=lambda: controller.open_page(
							  route='/options', source='/metainf', data=data))\
					.classes('w-2/3 self-center h-10 font-mono mt-6')

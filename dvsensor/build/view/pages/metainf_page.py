from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if not view.model.has_record_data:
		view.show_error('page not available')
		return
	back_button('input', view, on_click=lambda: view.model.clear())

	with ui.row().classes('w-full justify-center mt-24'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/3 h-auto '):
			with ui.column().classes('w-full place-content-center'):

				ui.label('Please check if the sequence record information is correct')\
					.classes('text-center text-lg font-mono')

				with ui.grid(columns=2).classes('place-self-start w-2/3'):
					ui.label('Name')\
						.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
					name_input = ui.input().classes('text-center text-base font-mono font-semibold')
					name_input.value = view.model.record_name
					name_input.bind_value(view.model, 'record_name')

					ui.label('ID')\
						.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
					id_input = ui.input().classes('text-center text-base font-mono font-semibold')
					id_input.value = view.model.record_id
					id_input.bind_value(view.model, 'record_id')

					ui.label('Description')\
						.classes('justify-self-end text-base font-mono font-semibold place-self-center mr-5')
					description_input = ui.input().classes('text-center text-base font-mono font-semibold')
					description_input.value = view.model.record_description
					description_input.bind_value(view.model, 'record_description')

				ui.button('ok',
						  on_click=lambda: view.open_page('options'))\
					.classes('w-2/3 self-center h-10 font-mono mt-6')

	footer()

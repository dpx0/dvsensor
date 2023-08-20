from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if view.model.sequence_record is None:
		view.show_error('page not available')
		return
	back_button('input', view, on_click=lambda: view.model.clear())

	with ui.row().classes('w-full justify-center mt-24'):
		with ui.card().classes('no-shadow border-[1px] rounded-xl p-8 w-2/3 h-64 '):
			with ui.column().classes('w-full place-content-center'):
				ui.label('Please check if the sequence record information is correct'
						 ).classes('text-center text-lg font-mono')


	footer()
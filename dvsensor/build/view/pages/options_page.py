from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()
	if not view.model.has_record_data:
		view.show_error('page not available')
		return
	back_button('metainf', view, on_click=lambda: view.model.reload_sequence_record())

	ui.label(view.model.record_name)
	ui.label(view.model.record_id)
	ui.label(view.model.record_description)
	ui.label(str(view.model.sequence))

	footer()

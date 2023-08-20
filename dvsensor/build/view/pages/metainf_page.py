from nicegui import ui
from ..base_elements import header, footer
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()

	if view.model.sequence_record is None:
		view.show_error('page not available')
		ui.button('go to main menu', on_click=lambda: view.open_page('start'))

	footer()
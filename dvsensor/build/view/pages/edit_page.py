from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	header()

	ui.label('Welcome to the Edit Page')
	ui.label(f'View reference {id(view)}')

	footer()
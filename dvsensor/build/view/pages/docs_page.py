from nicegui import ui
from ..base_elements import header, footer, back_button
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()
	back_button('start', view)

	ui.label('Welcome to the Help/Docs Page')
	ui.label(f'View reference {id(view)}')

	footer()


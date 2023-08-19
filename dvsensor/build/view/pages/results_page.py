from nicegui import ui
from ..base_elements import header, footer
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	set_colors()
	header()

	ui.label('Welcome to the Results Page')
	ui.label(f'View reference {id(view)}')

	footer()
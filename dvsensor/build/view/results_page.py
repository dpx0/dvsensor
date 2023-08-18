from nicegui import ui


def build(view, **kwargs) -> None:
	ui.label('Welcome to the Results Page')
	ui.label(f'View reference {id(view)}')

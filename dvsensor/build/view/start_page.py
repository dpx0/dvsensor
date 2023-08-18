from nicegui import ui


def build(view, **kwargs) -> None:
	ui.label('Welcome to the Start Page')
	ui.label(f'View reference {id(view)}')
	ui.button('Design Sensors',
			  on_click=view.controller.open_input_page
			  ).classes('w-32')


from nicegui import ui


def build(view, **kwargs) -> None:
	ui.label('Welcome to the Options Page')
	ui.label(f'View reference {id(view)}')
	ui.label(kwargs.get('sequence'))
	ui.button('Go',
			  on_click=lambda: view.controller.generate_sensors()).classes('w-32')
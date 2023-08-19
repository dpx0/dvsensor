from nicegui import ui
from .style import Colors


def header() -> None:
	with ui.header().style(f'background-color: {Colors.BACKGROUND}'
						   ).classes('items-center justify-between'):
		ui.label('DVSensor v1.0'
				 ).classes('text-2xl font-semibold tracking-wide font-mono'
				 ).style(f'color: {Colors.ACCENT}')


def footer() -> None:
	with ui.footer().style(f'background-color: {Colors.BACKGROUND}'):
		ui.label('FOOTER')


def back_button(previous_page: str, view) -> None:
	with ui.button('',
				   on_click=lambda: view.router.open(view.pages[previous_page]),
				   color=f'{Colors.ACCENT}'
				   ).props(f'no-caps flat round size=md'):
		ui.image("/assets/img-back01.png").classes('m-1')


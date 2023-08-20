from nicegui import ui
from .style import Colors
from typing import Callable, Optional


def header() -> None:
	with ui.header().style(f'background-color: {Colors.BACKGROUND}'
						   ).classes('items-center justify-between'):
		ui.label('DVSensor v1.0')\
			.classes('text-2xl font-semibold tracking-wide font-mono')\
			.style(f'color: {Colors.ACCENT}')


def footer() -> None:
	with ui.footer().style(f'background-color: {Colors.BACKGROUND}'):
		ui.label('FOOTER')


def back_button(previous_page: str, view, on_click: Optional[Callable] = None) -> None:
	def _on_click():
		if on_click:
			on_click()
		view.router.open(view.pages[previous_page])

	with ui.button('',
				   on_click=_on_click,
				   color=f'{Colors.ACCENT}')\
			.props(f'no-caps flat round size=md'):
		ui.image('/assets/img-back01.png').classes('m-1')


def show_dialog_box(message_lines: str | list[str], button_text: str = 'OK'):
	if type(message_lines) == str:
		message_lines = [message_lines]

	with ui.dialog() as dialog, ui.card():
		for msg in message_lines:
			ui.label(msg).classes('text-lg self-center font-mono')
		ui.button(button_text, on_click=dialog.close).classes('w-full self-center text-sm font-mono')
	dialog.open()

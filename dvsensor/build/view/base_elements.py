from nicegui import ui
from .style import Colors
from typing import Callable, Optional


def header() -> None:
	header = ui.header()\
		.style(f'background-color: {Colors.BACKGROUND}')\
		.classes('items-center justify-between')

	with header:
		ui.label('DVSensor v1.0')\
			.classes('text-2xl font-semibold tracking-wide font-mono')\
			.style(f'color: {Colors.ACCENT}')


def footer() -> None:

	with ui.footer().style(f'background-color: {Colors.BACKGROUND}'):
		...


def back_button(previous_page: str, view,
				on_click: Optional[Callable] = None, **kwargs) -> None:

	def _on_click():
		if on_click:
			on_click()
		view.router.open(view.pages[previous_page], **kwargs)

	with ui.button('', on_click=_on_click, color=f'{Colors.ACCENT}')\
			.props(f'no-caps flat round size=md'):
		ui.image('/assets/img-back01.png').classes('m-1')


def home_button(view):

	with ui.button('', on_click=lambda: view.router.open(view.pages['start']),
				   color=f'{Colors.ACCENT}')\
			.props(f'no-caps flat round size=md')\
			.classes('mt-4'):
		ui.image('/assets/img-home.png').classes('m-1')


def show_dialog_box(message_lines: str | list[str], button_texts: str | list[str] = 'OK'):
	if type(message_lines) == str:
		message_lines = [message_lines]

	if type(button_texts) == str:
		button_texts = [button_texts]

	with ui.dialog() as dialog, ui.card():
		for msg in message_lines:
			ui.label(msg).classes('text-lg self-center font-mono')
		with ui.row().classes('w-full'):
			for bt in button_texts:
				ui.button(bt, on_click=dialog.close).classes('w-full self-center text-sm font-mono')
	dialog.open()

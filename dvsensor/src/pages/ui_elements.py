from nicegui import ui
from typing import Any
from interface import Controller
from .style import Colors


def header() -> None:
	_header = ui.header()\
		.style(f'background-color: {Colors.BACKGROUND}')\
		.classes('items-center justify-between')

	with _header:
		ui.label('DVSensor v1.0')\
			.classes('text-2xl font-semibold tracking-wide font-mono')\
			.style(f'color: {Colors.ACCENT}')


def footer() -> None:
	with ui.footer().style(f'background-color: {Colors.BACKGROUND}'):
		...  # TODO: implement


def back_button(controller: Controller, route: str, source: str | None = None,
				data: dict[str, Any] | None = None) -> None:

	with ui.button('', on_click=lambda: controller.open_page(route=route, source=source, data=data),
				   color=f'{Colors.ACCENT}')\
			.props(f'no-caps flat round size=md'):
		ui.image('/assets/img-back01.png').classes('m-1')


def show_error_page(message: str, controller: Controller):
	ui.label(f"error: {message}").style(f'color: {Colors.RED}')
	ui.button('go to main menu', on_click=lambda: controller.open_page('/'))


def show_page_not_available(controller: Controller):
	show_error_page("page not available", controller)


def show_dialog_box(message_lines: str | list[str],
					button_texts: str | list[str] = 'OK') -> dict[str, ui.button]:
	# TODO: improve this function -> text should automatically be formated/placed/spaced correctly
	if type(message_lines) == str:
		message_lines = [message_lines]

	if type(button_texts) == str:
		button_texts = [button_texts]

	buttons = {}
	with ui.dialog() as dialog, ui.card():
		for msg in message_lines:
			ui.label(msg).classes('text-lg self-center font-mono')
		with ui.row().classes('w-full'):
			for bt in button_texts:
				button = ui.button(bt, on_click=dialog.close)\
					.classes('w-full self-center text-sm font-mono')
				buttons[bt] = button
	dialog.open()

	return buttons

from nicegui import ui
from functools import partial
from typing import Any
from interface import Controller
from .style import Colors
from .pagebuilder import page_builder


@page_builder()
def build(controller: Controller, data: dict[str, Any] | None) -> None:
	open_page = partial(controller.open_page, source='/')

	with ui.row().classes('w-full justify-center mt-32'):

		with ui.button('',
					   on_click=lambda: open_page('/input'),
					   color=f'{Colors.FOREGROUND}')\
				.props(f'no-caps text-color=primary')\
				.classes("rounded-xl p-8 w-1/4 h-64"):
			with ui.column():
				ui.label('Generate sensors for an mRNA transcript').classes('text-xl font-mono')
				ui.image("/assets/img-rna01.png").classes('w-1/4 h-1/4 mt-6 self-center')

		# with ui.button('', on_click=lambda: ..., color=f'{Colors.FOREGROUND}')\
		# 		.props(f'no-caps text-color=primary')\
		# 		.classes("rounded-xl p-8 w-1/4 h-64"):
		# 	with ui.column():
		# 		ui.label('Compare two mRNA transcripts').classes('text-xl font-mono')
		# 		ui.image("/assets/img-rna01db.png").classes('w-1/3 h-1/3 mt-6 self-center')

		with ui.button('',
					   on_click=lambda: open_page('/documentation'),
					   color=f'{Colors.FOREGROUND}')\
				.props(f'no-caps text-color=primary')\
				.classes("rounded-xl p-8 w-1/4 h-64"):
			with ui.column():
				ui.label('Help and Documentation').classes('text-xl font-mono')
				ui.image("/assets/img-book01.png").classes('w-1/3 h-1/3 mt-6 self-center')

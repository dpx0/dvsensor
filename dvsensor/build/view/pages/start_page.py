from nicegui import ui
from ..base_elements import header, footer
from ..style import Colors, set_colors


def build(view, **kwargs) -> None:
	view.controller.cancel_task()
	set_colors()
	header()

	with ui.row().classes('w-full justify-center mt-32'):

		with ui.button('',
					   on_click=lambda: view.router.open(view.pages['input']),
					   color=f'{Colors.FOREGROUND}')\
				.props(f'no-caps text-color=primary')\
				.classes("rounded-xl p-8 w-1/4 h-64"):
			with ui.column():
				ui.label('Generate sensors for an mRNA transcript').classes('text-xl font-mono')
				ui.image("/assets/img-rna01.png").classes('w-1/4 h-1/4 mt-6 self-center')

		with ui.button('', on_click=lambda: ..., color=f'{Colors.FOREGROUND}')\
				.props(f'no-caps text-color=primary')\
				.classes("rounded-xl p-8 w-1/4 h-64"):
			with ui.column():
				ui.label('Compare two mRNA transcripts').classes('text-xl font-mono')
				ui.image("/assets/img-rna01db.png").classes('w-1/3 h-1/3 mt-6 self-center')

		with ui.button('',
					   on_click=lambda: view.router.open(view.pages['documentation']),
					   color=f'{Colors.FOREGROUND}')\
				.props(f'no-caps text-color=primary')\
				.classes("rounded-xl p-8 w-1/4 h-64"):
			with ui.column():
				ui.label('Help and Documentation').classes('text-xl font-mono')
				ui.image("/assets/img-book01.png").classes('w-1/3 h-1/3 mt-6 self-center')

	footer()

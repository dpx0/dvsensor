"""
The page-routing functionality in this code was adapted from the nicegui example 'Single Page App'
by falkoschindler, see https://github.com/zauberzeug/nicegui/tree/main/examples/single_page_app
"""

from nicegui import ui
from .router import Router
from .pages.start_page import build as _build_start_page
from .pages.input_page import build as _build_input_page
from .pages.options_page import build as _build_options_page
from .pages.results_page import build as _build_results_page
from .pages.docs_page import build as _build_docs_page
from .pages.metainf_page import build as _build_metainf_page
from .style import Colors


class WebView:

	def __init__(self, model) -> None:
		self.model = model
		self.controller = None
		self.router = Router()
		self.pages = {}

	def open_page(self, page: str, **kwargs) -> None:
		page_builder = self.pages.get(page)

		if page_builder:
			self.router.open(page_builder, **kwargs)
		else:
			self.show_error(f"page {page} does not exist")

	def page_allowed(self, page_kwargs: dict) -> bool:
		if not page_kwargs.get('task_id') or not self.controller.task_controller.current_task_id or \
				page_kwargs.get('task_id') != self.controller.task_controller.current_task_id:
			return False
		return True

	def show_error(self, message: str) -> None:
		ui.label(f"error: {message}").style(f'color: {Colors.RED}')
		ui.button('go to main menu', on_click=lambda: self.open_page('start'))

	def build_pages(self) -> None:

		@ui.page('/')
		@ui.page('/{_:path}')
		def _build_pages() -> None:

			@self.router.add('/')
			def build_start_page(**kwargs) -> None:
				_build_start_page(self, **kwargs)
			self.pages['start'] = build_start_page

			@self.router.add('/documentation')
			def build_docs_page(**kwargs) -> None:
				_build_docs_page(self, **kwargs)
			self.pages['documentation'] = build_docs_page

			@self.router.add('/input')
			def build_input_page(**kwargs) -> None:
				_build_input_page(self, **kwargs)
			self.pages['input'] = build_input_page

			@self.router.add('/metainf')
			def build_metainf_page(**kwargs) -> None:
				_build_metainf_page(self, **kwargs)
			self.pages['metainf'] = build_metainf_page

			@self.router.add('/options')
			def build_options_page(**kwargs) -> None:
				_build_options_page(self, **kwargs)
			self.pages['options'] = build_options_page

			@self.router.add('/results')
			def build_results_page(**kwargs) -> None:
				_build_results_page(self, **kwargs)
			self.pages['results'] = build_results_page

			self.router.frame().classes('w-full')

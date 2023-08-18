"""
The page-routing functionality in this code was adapted from the nicegui example 'Single Page App'
by falkoschindler, see https://github.com/zauberzeug/nicegui/tree/main/examples/single_page_app
"""

from nicegui import ui
from .router import Router
from .start_page import build as _build_start_page
from .input_page import build as _build_input_page
from .options_page import build as _build_options_page
from .results_page import build as _build_results_page
from .edit_page import build as _build_edit_page


class WebView:

	def __init__(self, model) -> None:
		self.model = model
		self.controller = None
		self.router = Router()
		self.pages = {}

	def open_page(self, page: str, **kwargs) -> None:
		page_builder = self.pages.get(page)
		self.router.content.clear()
		if page_builder:
			page_builder(**kwargs)
		else:
			self.show_error(f"page {page} does not exist")

	@staticmethod
	def show_error(message: str) -> None:
		ui.label(f"error: {message}")

	def build_pages(self) -> None:

		@ui.page('/')
		@ui.page('/{_:path}')
		def _build_pages() -> None:

			@self.router.add('/')
			def build_start_page(**kwargs) -> None:
				_build_start_page(self, **kwargs)
			self.pages['start'] = build_start_page

			@self.router.add('/input')
			def build_input_page(**kwargs) -> None:
				_build_input_page(self, **kwargs)
			self.pages['input'] = build_input_page

			@self.router.add('/options')
			def build_options_page(**kwargs) -> None:
				_build_options_page(self, **kwargs)
			self.pages['options'] = build_options_page

			@self.router.add('/results')
			def build_results_page(**kwargs) -> None:
				_build_results_page(self, **kwargs)
			self.pages['results'] = build_results_page

			@self.router.add('/edit')
			def build_edit_page(**kwargs) -> None:
				_build_edit_page(self, **kwargs)
			self.pages['edit'] = build_edit_page

			self.router.frame().classes('w-full p-4 bg-gray-100')

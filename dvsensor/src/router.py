"""
This code was taken from the nicegui example 'Single Page App' by falkoschindler and
modified (changes indicated by * ). Minor changes regarding type annotations are not indicated.
See https://github.com/zauberzeug/nicegui/tree/main/examples/single_page_app for the original code.
"""
from typing import Awaitable, Callable, Dict, Any
from nicegui import background_tasks, ui
from interface import Controller, PageBuilderFn


class RouterFrame(ui.element, component='router_frame.js'):
	pass


class PageRouter:

	def __init__(self, page_builders: dict[str, PageBuilderFn]) -> None:
		self.routes: Dict[str, Callable] = {}
		self.content: ui.element | None = None
		self.register_pages(page_builders)  # * method call
		self.build_single_page_frame()  # * method call
		self.controller: Controller | None = None  # * added controller attribute

	# * removed add method

	# * added register_pages method
	def register_pages(self, pages: dict[str, PageBuilderFn]) -> None:
		for route, page_builder_fn in pages.items():
			self.routes[route] = page_builder_fn

	def open(self, path: str, source: str | None = None,
			 data: dict[str, Any] | None = None) -> None:
		# * replaced 'target' parameter with 'path', which can only be a string
		# * added data parameter
		self.controller.on_page_open(data)  # * method call

		# * always redirect to start page if requested page does not exist
		builder: PageBuilderFn = self.routes.get(path, self.routes['/'])

		async def build() -> None:
			with self.content:
				# * added JS code to remove all non-nicegui content (element class: 'non-ng')
				await ui.run_javascript(f'''
					var toRemove = document.getElementsByClassName("non-ng");
					Array.from(toRemove).forEach((element, i) => element.parentNode.removeChild(element));
					if (window.location.pathname !== "{path}") {{
						history.pushState({{page: "{path}"}}, "", "{path}");
					}}
				''', respond=False)
				# * builder functions are called with controller, source and data as arguments
				result = builder(self.controller, source, data)
				if isinstance(result, Awaitable):
					await result
		self.content.clear()
		background_tasks.create(build())

	def frame(self) -> ui.element:
		self.content = RouterFrame().on('open', lambda e: self.open(e.args))
		return self.content

	# * added build_single_page_frame method
	def build_single_page_frame(self) -> None:

		@ui.page('/')
		@ui.page('/{_:path}')
		def _build():
			self.frame().classes('w-full')

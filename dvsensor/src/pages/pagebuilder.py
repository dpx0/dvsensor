from typing import Callable, Any
from interface import Controller, PageBuilderFn
from .style import set_colors
from .ui_elements import header, footer, show_page_not_available


def page_builder(allow_from: str | list[str] | None = None) -> Callable[[Callable], PageBuilderFn]:

	def _decorator(build_fn: Callable[[Controller,  dict[str, Any] | None], None]) -> PageBuilderFn:
		def _build(controller: Controller, source: str | None, data: dict[str, Any] | None) -> None:
			set_colors()
			header()
			if isinstance(allow_from, str):
				_allow_from = [allow_from]
			else:
				_allow_from = allow_from

			if _allow_from is None or source in _allow_from:
				build_fn(controller, data)
			else:
				show_page_not_available(controller)
			footer()
		return _build

	return _decorator

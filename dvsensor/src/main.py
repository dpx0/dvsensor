import os
import sys
import logging
import importlib
import nicegui as ng
from router import PageRouter
from signals import SignalHandler
from controller import AppController
from interface import PageBuilderFn


def import_page_builder(module_name: str) -> PageBuilderFn:
	try:
		page_module = importlib.import_module('pages.' + module_name)
		return page_module.build
	except ModuleNotFoundError:
		logging.error(f"module pages.{module_name} not found")
		sys.exit(1)


class Application:

	def __init__(self) -> None:
		self_path = os.path.abspath(os.path.dirname(__file__))
		ng.app.add_static_files('/assets', self_path + '/assets')
		ng.app.on_exception(self.on_exception)

		page_builders = {
			'/': import_page_builder('page_start'),
			'/documentation': import_page_builder('page_docs'),
			'/input': import_page_builder('page_input'),
			'/metainf': import_page_builder('page_metainf'),
			'/options': import_page_builder('page_options'),
			'/results': import_page_builder('page_results'),
		}
		self.page_router = PageRouter(page_builders)
		self.controller = AppController(self.page_router, SignalHandler())

	@staticmethod
	def run() -> None:
		ng.ui.run(show=False, title='DVSensor v1.0', reload=False)

	@staticmethod
	def on_exception() -> None:
		ng.app.shutdown()


if __name__ in {'__main__', '__mp_main__'}:
	logging.basicConfig(level=logging.DEBUG,
	 					format='[%(levelname)s]\t %(funcName)s :: %(message)s')
	app = Application()
	app.run()

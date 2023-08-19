import nicegui as ng
from view import WebView
from model.datahandler import RNADataHandler
from controller import Controller
import os


class App:

	def __init__(self) -> None:
		self.model = RNADataHandler()
		self.view = WebView(self.model)
		self.controller = Controller(self.model, self.view)
		self.view.controller = self.controller

		self.root_path = os.path.abspath(os.path.dirname(__file__))
		ng.app.add_static_files('/assets', self.root_path + '/assets')

		ng.app.on_startup(self.on_startup)
		ng.app.on_shutdown(self.on_shutdown)
		ng.app.on_connect(self.on_connect)
		ng.app.on_disconnect(self.on_disconnect)
		ng.app.on_exception(self.on_exception)

	@staticmethod
	def on_startup() -> None:
		...

	@staticmethod
	def on_shutdown() -> None:
		...

	@staticmethod
	def on_connect() -> None:
		...

	@staticmethod
	def on_disconnect() -> None:
		...

	@staticmethod
	def on_exception() -> None:
		ng.app.shutdown()

	def start(self) -> None:
		self.view.build_pages()
		#ng.ui.run(show=False, reload=False)
		ng.ui.run(show=False)


def main() -> None:
	app = App()
	app.start()


if __name__ in {'__main__', '__mp_main__'}:
	main()

import nicegui as ng
from view import WebView
from controller import Controller


class App:

	def __init__(self) -> None:
		self.model = None
		self.view = WebView(self.model)
		self.controller = Controller(self.model, self.view)
		self.view.controller = self.controller

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
		ng.app.on_startup(self.on_startup)
		ng.app.on_shutdown(self.on_shutdown)
		ng.app.on_connect(self.on_connect)
		ng.app.on_disconnect(self.on_disconnect)
		ng.app.on_exception(self.on_exception)
		self.view.build_pages()
		ng.ui.run(show=False, reload=False)


def main() -> None:
	app = App()
	app.start()


if __name__ in {"__main__", "__mp_main__"}:
	main()


class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view

	def open_input_page(self) -> None:
		self.view.open_page('input')

	def process_input_sequence(self, sequence) -> None:
		self.view.open_page('options', sequence=sequence)




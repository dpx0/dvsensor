
class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view

	def process_fasta_sequence(self, sequence: str) -> None:
		self.view.open_page('options', sequence=sequence)

	def generate_sensors(self) -> None:
		self.view.open_page('results')



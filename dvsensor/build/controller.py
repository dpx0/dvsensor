import utils


class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view

	def handle_fasta_seq_input(self, user_input: str) -> None:
		seq_record = utils.read_fasta_string(user_input)
		self.model.load_sequence_record(seq_record)
		self.view.open_page('options')

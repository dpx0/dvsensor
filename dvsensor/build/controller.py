import utils


class Controller:

	def __init__(self, model, view) -> None:
		self.model = model
		self.view = view

	def handle_seq_file_upload(self, upload) -> None:
		input_str: str = upload.content.read().decode('UTF-8')
		seq_record = utils.read_fasta_string(input_str)
		print(seq_record)
		self.view.open_page('options', sequence=input_str)

	def handle_manual_seq_input(self, input_str: str) -> None:
		seq_record = utils.read_fasta_string(input_str)
		print(seq_record)
		self.view.open_page('options', sequence=input_str)





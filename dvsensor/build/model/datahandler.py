

class RNADataHandler:

	def __init__(self):
		self.sequence_record = None

	def clear(self):
		self.sequence_record = None

	def load_sequence_record(self, seq_record):
		self.sequence_record = seq_record

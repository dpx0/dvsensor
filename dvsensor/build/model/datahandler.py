class RNADataHandler:

	def __init__(self) -> None:
		self.init_sequence_record = None
		self.sequence = None
		self.record_name: str = ""
		self.record_id: str = ""
		self.record_description: str = ""
		self.has_record_data: bool = False

	def clear(self) -> None:
		self.init_sequence_record = None
		self.sequence = None
		self.record_name = ""
		self.record_id = ""
		self.record_description = ""
		self.has_record_data = False

	@property
	def sequence_length(self):
		if self.sequence:
			return len(self.sequence)

	def reload_sequence_record(self) -> None:
		self.load_sequence_record(self.init_sequence_record)

	def load_sequence_record(self, seq_record) -> None:
		self.init_sequence_record = seq_record
		self.sequence = seq_record.seq
		self.record_name = seq_record.name
		self.record_id = seq_record.id
		self.record_description = seq_record.description
		self.has_record_data = True

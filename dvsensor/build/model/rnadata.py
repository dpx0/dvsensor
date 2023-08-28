from dataclasses import dataclass, field
from Bio.SeqIO import SeqRecord


@dataclass
class RNAData:
	original_sequence_record: SeqRecord
	record_name: str = field(init=False)
	record_id: str = field(init=False)
	record_description: str = field(init=False)
	sequence: str = field(init=False)
	sequence_length: int = field(init=False)

	def __post_init__(self):
		self.record_name = self.original_sequence_record.name
		self.record_id = self.original_sequence_record.id
		self.record_description = self.original_sequence_record.description
		self.sequence = str(self.original_sequence_record.seq.transcribe())  # assume cDNA as input
		self.sequence_length = len(self.sequence)



import io
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def extract_sequence_data(sequence_record: SeqRecord) -> dict[str, str]:
	return {
		'sequence': str(sequence_record.seq.transcribe()),  # transcribe sequence (DNA or RNA) to RNA
		'name': sequence_record.name,
		'id': sequence_record.id
	}


def read_sequence(sequence: str) -> dict[str, str]:
	for fmt in ('fasta', 'genbank'):
		try:
			sequence_record: SeqRecord = SeqIO.read(io.StringIO(sequence), fmt)
			sequence_data: dict[str, str] = extract_sequence_data(sequence_record)
			if set(sequence_data['sequence']) - {'A', 'C', 'G', 'U'}:
				raise ValueError  # sequence contains letters other than A,C,G,U (after transcribing to RNA)
			return sequence_data
		except ValueError:
			pass  # reading sequence as single-entry FASTA failed, try GenBank format

	raise ValueError  # sequence is neither a (single-entry) FASTA nor GenBank record


import io
from Bio import SeqIO


def read_fasta_string(string: str):
	handle = io.StringIO(string)
	try:
		return SeqIO.read(handle, 'fasta')
	except ValueError:
		return None

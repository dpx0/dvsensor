import io
from Bio import SeqIO


def read_fasta_string(string: str):
	return SeqIO.read(io.StringIO(string), 'fasta')

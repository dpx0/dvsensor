import sys
import csv


def write_csv(output_file, header, rows):

	with open(output_file, "w") as fp:
		csvwriter = csv.writer(fp, dialect='excel', delimiter=',')
		csvwriter.writerow(header)
		for row in rows:
			csvwriter.writerow(row)


def read_csv(input_file):

	pass


def log_err(msg):
	print('[error] ', msg)
	sys.exit(1)


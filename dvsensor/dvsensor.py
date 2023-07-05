#!/usr/bin/env python3

"""
Authors:
	Daniel Prib      		<bytes@mailbox.org>
	Maximilian Leo Huber    <...>
	Saint Fischer 	 		<...>
Version: 1.0
Python Version: 3.11.3
Dependencies: biopython (1.81), numpy (1.25.0)
License: MIT License
"""

import os
import sys
from argparse import ArgumentParser
import lib.locate as locate

TRIPLETS = ({"CCA", "GCA", "UCA", "CAA", "CUA", "ACA"})
REGIONS = ({"5UTR", "CDS", "3UTR"})


def main(args):
	if not os.path.exists(args.input) or not os.path.isfile(args.input):
		log_err(f"input: file {args.input} not found")

	## if os.path.exists(args.output) and os.path.isfile(args.output):
		## log_err(f"output: file {args.output} already exists")

	if args.mode in ("locate", "all"):
		exclude = set(args.exclude)
		for tr in exclude:
			if tr not in TRIPLETS:
				log_err(f"exclude: unknown triplet {tr} " +
						"(options: CCA, GCA, UCA, CAA, CUA, ACA)")

		triplets = TRIPLETS - exclude
		if not triplets:
			log_err("exclude: at least one triplet is required")

		if not args.regions:
			regions = REGIONS
		else:
			regions = set(args.regions)
			for rg in regions:
				if rg not in REGIONS:
					log_err(f"[error]  regions: unknown region {rg} " +
							"(options: 5UTR, CDS, 3UTR)")

	if args.mode in ("generate", "all"):
		if args.ntleft <= 0:
			log_err("ntleft: must be greater than 0")
		if args.ntright <= 0:
			log_err("ntright: must be greater than 0")
	# TODO: process --blast parameter

	if args.mode == "locate":
		try:
			locateResults = locate.locate_triplets(args.input, triplets, regions)

		except (IOError, OSError) as err:
			log_err("could not read file '" + args.input + "': " + err.strerror)

		try:
			locate.write_output(args.output, locateResults)

		except (IOError, OSError) as err:
			log_err("could not write to file '" + args.output + "': " + err.strerror)

	elif args.mode == "generate":
		# TODO: implement 'generate' mode
		print("<debug> ntleft: ", args.ntleft)
		print("<debug> ntright: ", args.ntright)
	else:
		# TODO: implement 'all' mode
		print("<debug> triplets: ", triplets)
		print("<debug> regions : ", regions)
		print("<debug> ntleft: ", args.ntleft)
		print("<debug> ntright: ", args.ntright)


def log_err(msg):
	print('[error] ', msg)
	sys.exit(1)


if __name__ == "__main__":
	parser = ArgumentParser(description="")
	parser.add_argument("mode", choices=["locate", "generate", "all"])
	parser.add_argument("input", help="path to cDNA input file in FASTA format")
	parser.add_argument("output", help="path to output file in CSV format")
	locateArgs = parser.add_argument_group("locate")
	locateArgs.add_argument("--exclude", type=str, nargs="+", default=[],
							action="extend", metavar="triplets",
							help="list of triplets to be excluded (default = none) " +
								 "(options: CCA, GCA, UCA, CAA, CUA, ACA)")
	locateArgs.add_argument("--regions", type=str, nargs="+",
							action="extend", metavar="regions",
							help="list of transcript regions to be included " +
								 "(default = all) (options: 5UTR, CDS, 3UTR)")
	locateArgs.add_argument("--nosigp", action="store_true",
							help="disable prediction of signal peptides")
	generateArgs = parser.add_argument_group("generate")
	generateArgs.add_argument("--ntleft", metavar="nt", type=int, default=48,
							  help="number of nucleotides in sensor sequence " +
								   "left of stopcodon (default: %(default)s)")
	generateArgs.add_argument("--ntright", metavar="nt", type=int, default=48,
							  help="number of nucleotides in sensor sequence " +
								   "right of stopcodon (default: %(default)s)")
	generateArgs.add_argument("--blast", choices=["ncbi", "local"],
							  help="")
	# TODO: add help text for --blast parameter
	# TODO: add additional parameters for BLAST-querying (e.g. path to local database)
	main(parser.parse_args())


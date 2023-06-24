#!/usr/bin/env python3

"""
Authors:
	Daniel Prib      <bytes@mailbox.org>
	M. H.            <...>
	Saint Fischer 	 <...>
Version: 1.0
Python Version: 3.11.3
Dependencies: biopython (1.81), numpy (1.25.0)
License: MIT License
"""

import os
import sys
from argparse import ArgumentParser


TRIPLETS = ({"CCA", "GCA", "UCA", "CAA", "CUA", "ACA"})
REGIONS = ({"5UTR", "CDS", "3UTR"})


def main(args):
	
	if not os.path.exists(args.input) or not os.path.isfile(args.input):
		print(f"input: file {args.input} not found")
		sys.exit(1)
	
	if os.path.exists(args.output) and os.path.isfile(args.output):
		print(f"output: file {args.output} already exists")
		sys.exit(1)
	
	if args.mode in ("locate", "all"):
		exclude = set(args.exclude)
		for tr in exclude:
			if tr not in TRIPLETS:
				print(f"exclude: unknown triplet {tr} " + 
				       "(options: CCA, GCA, UCA, CAA, CUA, ACA)")
				sys.exit(1)
		triplets = TRIPLETS - exclude
		if not triplets:
			print("exclude: at least one triplet is required")
			sys.exit(1)
		
		if not args.regions:
			regions = REGIONS
		else:
			regions = set(args.regions)
			for rg in regions:
				if rg not in REGIONS:
					print(f"regions: unknown region {rg} " +
					       "(options: 5UTR, CDS, 3UTR)")
					sys.exit(1)
				
	if args.mode in ("generate", "all"):
		if args.ntleft <= 0:
			print("ntleft: must be greater than 0")
			sys.exit(1)
		if args.ntright <= 0:
			print("ntright: must be greater than 0 ")
			sys.exit(1)
		# TODO: process --blast parameter
			
	if args.mode == "locate":
		# TODO: implement 'locate' mode
		print("<debug> triplets: ", triplets)
		print("<debug> regions : ", regions)
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
		

if __name__ == "__main__":
	parser = ArgumentParser(description="")
	parser.add_argument("mode", choices= ["locate", "generate", "all"])
	parser.add_argument("input",
						help="path to cDNA input file in FASTA format")
	parser.add_argument("output",
						help="path to output file in CSV format")
	locateArgs = parser.add_argument_group("locate")
	locateArgs.add_argument("--exclude", type=str, nargs="+", default=[],
							action="extend", metavar="triplets",
							help="list of triplets to be excluded (default = none) " +
						         "(options: CCA, GCA, UCA, CAA, CUA, ACA)")
	locateArgs.add_argument("--regions", type=str, nargs="+", 
							action="extend", metavar="regions",
							help="list of transcript regions to be included "+
								 "(default = all) (options: 5UTR, CDS, 3UTR)")
	locateArgs.add_argument("--nosigp", action="store_true", 
							help="disable prediction of signal peptides")
	generateArgs = parser.add_argument_group("generate")
	generateArgs.add_argument("--ntleft", metavar="nt", type=int, default=48,
							help="number of nucleotides in sensor sequence "+
							     "left of stopcodon (default: %(default)s)")
	generateArgs.add_argument("--ntright", metavar="nt", type=int, default=48,
							help="number of nucleotides in sensor sequence "+
							     "right of stopcodon (default: %(default)s)")
	generateArgs.add_argument("--blast", choices=["ncbi", "local"],
							help="") 
	# TODO: add help text for --blast parameter
	# TODO: add additional parameters for BLAST-querying (e.g. path to local database)
	main(parser.parse_args())

#!/usr/bin/env python

import optparse


def main():
	#Parse Command Line
	parser = optparse.OptionParser(description='A small Python3 script to get Krona outputs from metaPhLan results',

	epilog="Adapted from Author: Bertalan Takács, @TakacsBertalan on GitHub")
	parser.add_option( '-p', '--profile', dest='profile', default='', action='store', help='Name of the metaPhLan3 taxonomic output file. Alternatively you can specify a folder path and a search string' )
	parser.add_option( '-k', '--krona', dest='krona', default='krona.out', action='store', help='the Krona output file name [krona.out]' )
	
	( options, spillover ) = parser.parse_args()
	if options.profile == "":
		print("Please specify an input metaplan profile file")
		parser.print_help()
		return
	if options.profile != "":
		convert_metaphlan_output(options.profile, options.krona)	



def convert_metaphlan_output(input_file, output_file):
	metaphlan_result = {}
	with open(input_file, "r") as inp:
		for line in inp:
			if line[0] != "#":
				line_list = line.rstrip().split("\t")
				if "s__" in line_list[0] or "UNCLASSIFIED" in line_list[0]:
					metaphlan_result[line_list[0]] = line_list[2]

	with open(output_file, "w") as outp:
		for key in metaphlan_result:
			taxonomy = key.split("|")
			taxonomy = [x[3:] if x != "UNCLASSIFIED" else x for x in taxonomy]
			outp.write(metaphlan_result[key] + "\t" + "\t".join([x.replace("_", " ", 2).replace("_","-") for x in taxonomy]) + "\n")


main()

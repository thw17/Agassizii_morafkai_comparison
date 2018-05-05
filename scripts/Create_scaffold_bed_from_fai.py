import argparse
import csv


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--fai", required=True,
		help="Input .fai file.")

	parser.add_argument(
		"--outfile", required=True,
		help="Output file (with full path).")

	parser.add_argument(
		"--scaffold_name", required=True,
		help="Name of scaffold in fai file")

	args = parser.parse_args()
	return args


def main():
	args = parse_args()

	with open(args.fai, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		for i in reader:
			if i[0] == args.scaffold_name:
				temp_line = [i[0], 0, i[1]]
				break

	with open(args.outfile, "w") as o:
		o.write(
			"{}\t{}\t{}".format(temp_line[0], temp_line[1], temp_line[2]))


if __name__ == "__main__":
	main()

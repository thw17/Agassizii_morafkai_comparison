import argparse
import csv


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--fai", required=True,
		help="Input .fai file.")

	parser.add_argument(
		"--out_prefix", required=True,
		help="Output prefix for bed files.")

	parser.add_argument(
		"--chunks", type=int, default=10,
		help="Number of approximately equal sized chunks to "
		"break genome into.")

	args = parser.parse_args()
	return args


def main():
	args = parse_args()

	total_length = 0

	with open(args.fai, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		for i in reader:
			total_length += int(i[1])

	chunk_size = float(total_length) / args.chunks

	max_chunk = chunk_size

	chunk_dict = {}
	for i in range(1, args.chunks + 1):
		chunk_dict[i] = []

	with open(args.fai, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		chunk_number = 1
		running_len = 0
		for i in reader:
			chunk_dict[chunk_number].append([i[0], 0, i[1]])
			running_len += int(i[1])
			if running_len >= max_chunk:
				max_chunk += chunk_size
				chunk_number += 1

	for i in chunk_dict:
		out_name = "{}_chunk{}.bed".format(
			args.out_prefix, i)
		with open(out_name, "w") as o:
			writer = csv.writer(o, delimiter="\t")
			for k in chunk_dict[i]:
				writer.writerow(k)


if __name__ == "__main__":
	main()

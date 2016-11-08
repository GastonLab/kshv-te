import csv
import sys
import argparse

from collections import defaultdict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', help="Output file name")
    parser.add_argument('-i', '--inputs', help='Input file name with list of Salmon quant files and paths')
    args = parser.parse_args()

    transcript_counts = defaultdict(dict)
    samples = list()

    with open(args.input, 'r') as inputs:
        for line in inputs:
            temp = line.split()
            sample = temp[0]
            infile = temp[1]
            samples.append(sample)

            with open(infile, 'r') as salmonfh:
                reader = csv.reader(salmonfh, dialect='excel-tab')
                header = reader.next()
                for row in reader:
                    transcript_counts[row[0]][sample] = int(row[4])

    with open(args.output, 'w') as outfile:
        outfile.write("Transcript")
        for sample in samples:
            outfile.write("\t{}".format(sample))
        outfile.write("\n")

        for transcript in transcript_counts.keys():
            outfile.write("{}".format(transcript))
            for sample in samples:
                outfile.write("\t{}".format(transcript_counts[transcript][sample]))
            outfile.write("\n")

import argparse
import sys
import collections

def main():
    args = get_args()
    genes_file = args.genes
    exon_file = args.exones
    output = args.output
    # Create set of all genes in the 273 genes (they are 250 only).
    genes = set()
    with open(genes_file, 'r') as data_in:
        for line in data_in:
            genes.add(line.split()[0])
    with open(exon_file, 'r') as data_in, open(output, 'w') as data_out:
        for line in data_in:
            gene = line.split()[5]
            if gene in genes:
                data_out.write(line)



def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for more info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Select fields from GTF file',
                                     add_help=True, )

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
    parser.add_argument('-g', "--genes",  help="<Required> genes file", metavar="FOO.tsv", required=True)
    parser.add_argument('-o', "--output", help="<Required> output exones file", metavar="OUT.tsv", required=True)
    parser.add_argument('-e', '--exones', help="<Required> exones file", required=True, dest='exones')


    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()

import argparse
import sys
import collections

def main():
    args = get_args()
    input_file = args.input
    features = args.features
    output = args.output
    feature_dictionary = {el:"." for el in features}
    first = True
    column_names = "chr\tstart\tend\t{}\n".format("\t".join(features))
    with open(input_file, 'r') as data_in, open(output, 'w') as data_out:
        line_split = []
        for line in data_in:
            if not line.startswith("#"):
                line_split = line.split("\t")
                if line_split[2] == "exon":
                    # print(line_split)
                    for feature in line_split[8].split(";"):
                        feature = feature.strip()
                        feature_split = feature.split(" ")
                        if feature_split[0] in features:
                            feature_dictionary[feature_split[0]] = feature_split[1].replace('"', '')
                    sorted_dict = collections.OrderedDict(feature_dictionary)
                    sorted_dict_values = "\t".join(sorted_dict.values())
                    if first:
                        data_out.write(column_names)
                        first = False
                        data_out.write(f'{line_split[0]}\t{line_split[3]}\t{line_split[4]}\t{sorted_dict_values}\n')
                    else:
                        data_out.write(f'{line_split[0]}\t{line_split[3]}\t{line_split[4]}\t{sorted_dict_values}\n')


# gene_id transcript_id gene_name exon_id

def gen_value(my_feature):
    if my_feature:
        feature_split = my_feature.split(" ")
        return feature_split[1].replace('"', '')


def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for more info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Select fields from GTF file',
                                     add_help=True, )

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
    parser.add_argument('-i', "--input",  help="<Required> GTF file", metavar="FOO.gtf", required=True)
    parser.add_argument('-o', "--output", help="<Required> selectd features from gtf file", metavar="OUT.tsv", required=True)
    parser.add_argument('-f', '--features', nargs='+', required=True, dest='features')


    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()

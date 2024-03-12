import sys


bench = sys.argv[1]
coverage = sys.argv[2]




bench_dic = {}
with open(bench, 'r') as data_in:
    for line in data_in:
        line_split = line.split()
        if len(line_split) < 6:
            line_len = len(line_split)
            missing_values = 6 - line_len + 1
            line_split.extend(['0']*missing_values)
        if 'NaN' in line_split:
            line_split = [ item if item != 'NaN' else '0' for item in line_split ]

        bench_dic[line_split[0]] = line_split[1:]
with open(coverage, 'r') as data_in:
    for line in data_in:
        gene = line.split()[4]
        try:
            values = "\t".join(bench_dic[gene])
        except Exception as e:
            print(gene, bench_dic[gene], e)

        print(f"{line.strip()}\t{values}")

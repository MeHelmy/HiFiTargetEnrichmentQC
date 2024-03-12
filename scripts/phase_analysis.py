#!/usr/bin/env python
import sys
from dataclasses import dataclass

@dataclass
class Gene:
    phase_blocks: set
    total_variant: int = 0
    total_het: int = 0
    phased_variant: int = 0


# use stdin if it's full
if not sys.stdin.isatty():
    input_stream = sys.stdin

# otherwise, read the given filename
else:
    try:
        input_filename = sys.argv[1]
    except IndexError:
        message = 'need filename as first argument if stdin is not full'
        raise IndexError(message)
    else:
        input_stream = open(input_filename, 'rU')

gene_stat = {}
for line in input_stream:
    # print(line)
    gene, info, gt = line.split()

    if gene in gene_stat:
        phased = 1 if "|" in gt else 0
        gene_stat[gene].total_variant = gene_stat[gene].total_variant + 1
        het = 1 if "0" in gt.split(':')[0] else 0
        gene_stat[gene].total_het = gene_stat[gene].total_het + het
        if phased:
            gene_stat[gene].phased_variant = gene_stat[gene].phased_variant + phased
            gene_stat[gene].phase_blocks.add(gt.split(':')[-1])
    else:
        phased = 1 if "|" in gt else 0
        het = 1 if "0" in gt.split(':')[0] else 0
        phase_block = set(gt.split(':')[-1]) if phased else set()
        gene_stat[gene] = Gene(phase_blocks=phase_block, total_variant=1, total_het=het, phased_variant=phased)
first=1
for k, v in gene_stat.items():
    if first:
        print("Gene\tPhase_blocks\tTotal_variants\tHetero_variants\tPhased_variants\tPercentage_phased")
        first=0
    try:
        phased_percentage =  0 if v.phased_variant == 0 else (v.phased_variant/v.total_het)*100 or 0
        print(f"{k}\t{len(v.phase_blocks)}\t{v.total_variant}\t{v.total_het}\t{v.phased_variant}\t{phased_percentage:.2f}")
    except Exception as e:
        print("wrong, can not divide by zero", file=sys.stderr)
        print(f"{k}\t{len(v.phase_blocks)}\t{v.total_variant}\t{v.total_het}\t{v.phased_variant}", file=sys.stderr)
        print(e, file=sys.stderr)
        exit(1)

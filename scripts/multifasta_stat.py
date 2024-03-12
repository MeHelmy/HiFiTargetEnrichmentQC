import sys, re
from collections import Counter
# import regex
import pyfastx

input = sys.argv[1]
homo_length = int(sys.argv[2])
bed = sys.argv[3]

# fa = pyfastx.Fastx(input)
# for name,seq,comment in fa:
#     print(f"Name:{name}\nSequence:{seq}\nComment:{comment}")
#     print(seq.gc_content)

bed_genes = {}

with open(bed, 'r') as data_in:
    for line in data_in:
        chr, start, end, gene, length = line.split()
        bed_genes[gene] = [chr, start, end, length]

print(f"Chr\tStart\tEnd\tLength\tGene\tA\tT\tC\tG\tpolymer{homo_length}+\tGC%\tAG%")
for name, seq in pyfastx.Fasta(input, build_index=False):
    nucleotide_count = Counter(seq)
    # polymer_length = sum([len(m.group()) for m in re.finditer(r'([ACGT])\1{9,}', seq) if len(m.group()) >= homo_length])
    polymer_length = sum([len(m.group()) for m in re.finditer(r"([ACGT])\1{"+str(homo_length-1)+",}", seq) if len(m.group()) >= homo_length])
    location = "\t".join(bed_genes[name])
    GC = (int(nucleotide_count['G']) + int(nucleotide_count['C']))/float(bed_genes[name][3])*100
    AG = (int(nucleotide_count['G']) + int(nucleotide_count['A']))/float(bed_genes[name][3])*100
    print(f"{location}\t{name}\t{nucleotide_count['A']}\t{nucleotide_count['T']}\t{nucleotide_count['C']}\t{nucleotide_count['G']}\t{polymer_length}\t{GC}\t{AG}")


# for seq in pyfastx.Fasta(input):
#     print(seq.name)
#     print(seq.seq)

# fa = pyfastx.Fasta(input)
# print(len(fa))
# print(fa.size)
# print(fa.gc_content)
# print(fa.composition)

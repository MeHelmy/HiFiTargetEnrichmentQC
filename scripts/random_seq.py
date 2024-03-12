#!/usr/bin/env python

import sys
import random
from Bio import SeqIO

sequence = sys.argv[1]
seq_length = int(sys.argv[2])
seq_type = sys.argv[3]


def main():
    with open(sequence) as handle:
        for record in SeqIO.parse(handle, seq_type):
            if len(record.seq) == seq_length:
                print(f">{record.id}\n{record.seq}")
            elif  len(record.seq) > seq_length:
                base = len(record.seq)//2
                new_seq = cut_read_recusrsive(record.seq, seq_length, base)
                print(f">{record.id}\n{new_seq}")

def cut_read_recusrsive(seq, read_leng, base):
    base_rand = random.randint(0, base)
    if (base_rand + read_leng) <= len(seq):
        myseq = seq[base_rand:(base_rand + read_leng)]
        return(myseq)
    else:
        return(cut_read_recusrsive(seq, read_leng, base))


if __name__ == "__main__":
    main()
            # base = len(record.seq)//2
            # base_rand = random.randint(0, base)
            # if base_rand + seq_length <= len(record.seq):
            #     print(f">{record.id}\n{record.seq[base_rand: base_rand + seq_length]}")

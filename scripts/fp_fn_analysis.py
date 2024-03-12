#!/usr/bin/env python
import sys
from dataclasses import dataclass




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
    # Peocess the bed file 

#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from pathlib import Path

def combine_sequences(raw_data_directory, output_file):
    p = Path(raw_data_directory)
    fasta_files = p.glob('*.fasta')
    for filename in fasta_files:
        for line in filename.open():
            if line.startswith('>'):
                pipe_position = line.index('|')
                slash_position = line.index('/')
                line = '>' + line[slash_position + 1:pipe_position]
            if line.endswith('\x13'):
                line = line[:-1]
            if not line.endswith('\n'):
                line += '\n'
            output_file.write(line)
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert GISAID sequences to Nextstrain input')
    parser.add_argument('raw_data_directory', help='Directory containing sequences from GISAID (with GISAID header format)')
    parser.add_argument('output_file', default='data/sequences.fasta', type=argparse.FileType('w'), help='Output filename (sequences.fasta for Nextstrain)')
    args = parser.parse_args()
    combine_sequences(args.raw_data_directory, args.output_file)

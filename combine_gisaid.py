#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from io import StringIO
from pathlib import Path
import re

from Bio import SeqIO

def combine_sequences(raw_data_directory, output_file):
    p = Path(raw_data_directory)
    fasta_files = p.glob('*.fasta')
    records = []
    keys_seen = set()
    for filename in fasta_files:
        buffer = StringIO()
        for line in filename.open():
            if line.startswith('>'):
                line = line.replace(' ', '') # make sure we get HongKong instead of Hong Kong
            if line.endswith('\x13'):
                line = line[:-1]
            if not line.endswith('\n'):
                line += '\n'
            buffer.write(line)
        buffer.seek(0)

        seq_record = SeqIO.read(buffer, 'fasta')
        if len(seq_record.seq) < 20000:
            # skip short sequences
            continue
        seq_record.id = re.sub(r'^[^/]*/([^|]*)\|.*', r'\1', seq_record.id)
        seq_record.id = seq_record.id.replace(' ', '')
        print(seq_record.id)
        if seq_record.id in keys_seen:
            continue
        keys_seen.add(seq_record.id)
        records.append(seq_record)
        # for line in filename.open():
        #     if line.startswith('>'):
        #         line = line.replace(' ', '')  # make sure we get HongKong instead of Hong Kong
        #         pipe_position = line.index('|')
        #         slash_position = line.index('/')
        #         line = '>' + line[slash_position + 1:pipe_position]

        #     output_file.write(line)
    SeqIO.write(records, output_file, 'fasta')
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert GISAID sequences to Nextstrain input')
    parser.add_argument('raw_data_directory', help='Directory containing sequences from GISAID (with GISAID header format)')
    parser.add_argument('output_file', nargs='?', default='data/sequences.fasta', type=argparse.FileType('w'), help='Output filename (sequences.fasta for Nextstrain)')
    args = parser.parse_args()
    combine_sequences(args.raw_data_directory, args.output_file)

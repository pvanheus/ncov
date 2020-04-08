#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from collections import Counter
from io import StringIO
from operator import itemgetter
from pathlib import Path
import re

from Bio import SeqIO,Seq
from Bio.Alphabet import generic_dna


def get_clean_file(filename):
    buffer = StringIO()
    for line in open(filename):
        if line.startswith('>'):
            line = line.replace(' ', '') # make sure we get HongKong instead of Hong Kong
        if line.endswith('\x13'):
            line = line[:-1]
        if not line.endswith('\n'):
            line += '\n'
        buffer.write(line)
    buffer.seek(0)
    return buffer


def combine_sequences(raw_data_location, output_file, min_length=20000, max_gaps_perc=5.0):
    p = Path(raw_data_location)
    records = []
    if p.is_dir():
        fasta_files = p.glob('*.fasta')
        buffer = StringIO()
        for filename in fasta_files:
            buffer.write(get_clean_file(filename))
    else:
        buffer = get_clean_file(raw_data_location)

    buffer.seek(0)
    keys_seen = set()
    for seq_record in SeqIO.parse(buffer, 'fasta'):
        if len(seq_record.seq) < min_length:
            # skip short sequences
            continue

        seq_record.id = re.sub(r'^[^/]*/([^|]*)\|.*', r'\1', seq_record.id)
        seq_record.id = seq_record.id.replace(' ', '')

        if seq_record.id in keys_seen:
            continue
        keys_seen.add(seq_record.id)
        seq_record.seq = Seq.Seq(str(seq_record.seq).upper().replace('-', 'N'), alphabet=generic_dna)
        gap_percentage = seq_record.seq.count('N') / len(seq_record.seq) * 100.0
        if gap_percentage > max_gaps_perc:
            continue

        base_counts = ' '.join([ f'{e[0]}:{e[1]}' for e in sorted(Counter(str(seq_record.seq)).items(), key=itemgetter(1), reverse=True)])
        print(seq_record.id, round(gap_percentage), base_counts, sep='\t')
        records.append(seq_record)

    SeqIO.write(records, output_file, 'fasta')
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert GISAID sequences to Nextstrain input')
    parser.add_argument('--min_seq_len', type=int, default=20000, help='Minimum length of the sequence')
    parser.add_argument('--max_gap_percentage', type=int, default=5.0, help='Maximum percentage of the sequence that may be gaps (- or N characters)')
    parser.add_argument('raw_data', help='File or Directory containing sequences from GISAID (with GISAID header format)')
    parser.add_argument('output_file', nargs='?', default='data/sequences.fasta', type=argparse.FileType('w'), help='Output filename (sequences.fasta for Nextstrain)')
    args = parser.parse_args()
    combine_sequences(args.raw_data, args.output_file, args.min_seq_len, args.max_gap_percentage)

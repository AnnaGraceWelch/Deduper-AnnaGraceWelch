#!/usr/bin/env python

import argparse
import re
from IPython import embed

def get_args():
    parser = argparse.ArgumentParser(description='program to remove PCR duplicates from a SAM file')
    parser.add_argument('-f', '--filename', help="sorted input SAM file remove duplicates from")
    parser.add_argument('-o', '--output', help='output file to write records minus duplicates')
    parser.add_argument('-d', '--dup_file', help='output file to write duplicate records to')
    parser.add_argument('-s', '--summary', help='output file to write summary of program results to')
    parser.add_argument('-u', '--umi_file', help='file containing list of known UMIs')
    return parser.parse_args()

def known_umis(fh) -> dict:
    '''Takes a file handle as input. Opens the file, appends each known UMI in the file to a dictionary. Returns set of UMIs.'''
    known = set()
    with open(fh, 'r') as input:
        for line in input:
            line = line.strip()
            known.add(line)
    return known

def parse_bitwise(flag: int) -> str:
    '''Takes a bitwise flag integer value as input. Checks whether sequence is reverse complement or not.
    If yes, returns the string 'minus'. If no, returns the string 'plus'.'''
    if ((flag & 16) == 16):
        strand = 'minus'
    else:
        strand = 'plus'
    return strand


def adjust_pos(cigar: str, pos: int, strand: str):
    '''This function takes a CIGAR string, position, and strandedness as input. It adjusts the 5' start position based on
   soft-clipping and strandedness. Returns the adjusted position.'''
    #find soft-clipping
    if strand == 'minus':
        m = re.findall(r'(\d+)[M|D|N]', cigar)
        for match in m:
            pos = pos + int(match)
        if 'S' in cigar:
            if cigar.endswith('S'):
                s = re.findall(r'(\d+)S', cigar)
                clipped = s[-1]
                #clipped = s.group(1)
                pos = pos + int(clipped)
    else:
        if 'S' in cigar:
            s = re.search(r'(\d+)S', cigar)
            clipped = s[1]
            pos = pos - int(clipped)
    return pos

args = get_args()

unknown_umis = 0
duplicates = 0
known = known_umis(args.umi_file)
records = 0

current_chrom = None
seen = {}

with open(args.filename, 'r') as input, open(args.output, 'w') as output, open(args.dup_file, 'w') as dupfile:
    for line in input:
        if line.startswith('@'):
            output.write(line)
        else:
            record = line.strip().split()
            chrom = record[2]
            if chrom != current_chrom: 
                current_chrom = chrom
                seen.clear()
            umi = line.split(':')[7]
            umi = umi.split()[0]
            if umi in known:
                pos = record[3]
                cigar = record[5]
                strand = parse_bitwise(int(record[1]))
                pos = adjust_pos(cigar, int(pos), strand)
                if (umi, pos, strand) in seen:
                    seen[(umi, pos, strand)] += 1
                    duplicates += 1
                    records += 1
                    dupfile.write(line)
                else:
                    output.write(line)
                    seen[(umi, pos, strand)] = 1
                    records += 1
            else: 
                unknown_umis += 1
                records += 1

with open(args.summary, 'w') as summary:
    summary.write(f'Number of duplicates removed: {duplicates}\n')
    summary.write(f'Percent of records that were duplicates: {(duplicates / records) * 100}%\n')
    summary.write(f'Number of unknown UMIs: {unknown_umis}\n')
    summary.write(f'Percent of records with unknown_umis: {(unknown_umis / records) * 100}%\n')
    summary.write(f'Number of records removed from original file: {unknown_umis + duplicates}\n')
    summary.write(f'Percent of records removed from original file: {((unknown_umis + duplicates) / records) * 100}%\n')
    summary.write(f'Percent of records retained: {100 -(((unknown_umis + duplicates) / records) * 100)}%')


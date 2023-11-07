#!/usr/bin/env python

import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description='program to remove PCR duplicates from a SAM file')
    parser.add_argument('-f', '--filename', help="sorted input SAM file remove duplicates from")
    parser.add_argument('-o', '--output', help='output file to write records minus duplicates')
    #parser.add_argument('-d', '--dup_file', help='output file to write duplicate records to')
    #parser.add_argument('-s', '--summary', help='output file to write summary of program results to')
    parser.add_argument('-u', '--umi_file', help='file containing list of known UMIs')
    return parser.parse_args()

def known_umis(fh) -> dict:
    '''Takes a file handle as input. Opens the file, appends each known UMI in the file to a dictionary. Returns set of UMIs.'''
    #create a set to hold known umis
    known = set()
    #open file with known umis
    with open(fh, 'r') as input:
        #loop through  file
        for line in input:
            #strip white space from each line
            line = line.strip()
            #add that umi to set 
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
    #if we're on the minus strand
    if strand == 'minus':
        #find all occurrences of matches, deletions, and skipped regions
        m = re.findall(r'(\d+)[M|D|N]', cigar)
        for match in m:
            # add each match to pos
            pos = pos + int(match)
        # check for soft-clipping
        if 'S' in cigar:
            #if soft-clipping right-handed
            if cigar.endswith('S'):
                # find all occurrences of soft-clipping:
                s = re.findall(r'(\d+)S', cigar)
                #capture only the last occurrence
                clipped = s[-1]
                #clipped = s.group(1)
                #add soft-clipping to pos
                pos = pos + int(clipped)
    else:
        # check for soft-clipping if we're on plus strand
        if 'S' in cigar:
            # search for left-hand soft-clipping
            if re.search(r'^(\d+)S', cigar):
                #capture value for soft-clipping
                s = re.search(r'^(\d+)S', cigar)
                clipped = s[1]
                #subtract soft-clipping from position
                pos = pos - int(clipped)
    return pos


# get input SAM file, name of output file, and file with known UMIs
args = get_args()

# initialize counters for unknown umis, duplicates, and overall records
unknown_umis = 0
duplicates = 0
records = 0
# populate set with known umis
known = known_umis(args.umi_file)

#initialize variable for the current chromosome
current_chrom = None
#initialize dictionary to hold tuples of seen (umis, positions, strandedness)
seen = {}

#open input file and output file
with open(args.filename, 'r') as input, open(args.output, 'w') as output, open('dupfile_sam', 'w') as dupfile:
    #loop through file
    for line in input:
        #if it's a header line
        if line.startswith('@'):
            #write header to output
            output.write(line)
        else:
            #split line on white space
            record = line.strip().split()
            #get chromosome field
            chrom = record[2]
            #if we've moved on to the next chromosome
            if chrom != current_chrom: 
                #set current chromosome to next chromosome
                current_chrom = chrom
                #empty dictionary of seen umis, pos, strand
                seen.clear()
            #split line on ':'
            umi = line.split(':')[7]
            #split to get umi value
            umi = umi.split()[0]
            #check if umi is known umi
            if umi in known:
                #get pos of record
                pos = record[3]
                #get cigar string of record
                cigar = record[5]
                #parse bitwise flag to get strandedness
                strand = parse_bitwise(int(record[1]))
                #adjust position 
                pos = adjust_pos(cigar, int(pos), strand)
                #if this is a duplicate
                if (umi, pos, strand) in seen:
                    #increment counter in dictionary
                    seen[(umi, pos, strand)] += 1
                    #increment counter of duplicates
                    duplicates += 1
                    #increment counter of total records
                    records += 1
                    #write duplicate to file of duplicates 
                    dupfile.write(line)
                else:
                    #it's not a duplicate so write record to output
                    output.write(line)
                    #add umi, pos, strand to dictionary
                    seen[(umi, pos, strand)] = 1
                    #increment counter for total records
                    records += 1
            else: 
                #umi is unknown
                #increment counter for unknown umis
                unknown_umis += 1
                #increment counter for total records
                records += 1

#open file to write summary to
with open('summary.txt', 'w') as summary:
    summary.write(f'Total Records: {records}\n')
    summary.write(f'Number of duplicates removed: {duplicates}\n')
    summary.write(f'Percent of records that were duplicates: {(duplicates / records) * 100}%\n')
    summary.write(f'Number of unknown UMIs: {unknown_umis}\n')
    summary.write(f'Percent of records with unknown_umis: {(unknown_umis / records) * 100}%\n')
    summary.write(f'Number of records removed from original file: {unknown_umis + duplicates}\n')
    summary.write(f'Percent of records removed from original file: {((unknown_umis + duplicates) / records) * 100}%\n')
    summary.write(f'Percent of records retained: {100 -(((unknown_umis + duplicates) / records) * 100)}%')


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:05:42 2021

@author: cfos
"""

from Bio import SeqIO
import os
import argparse
import sys
import multiprocessing as mp
import tqdm

#%%
def filter_alignment(record, gap_threshold, n_threshold, cumulative_threshold, outfile, write_filtered):
    flat_seq = str(record.seq)
    seq_name = record.id
    gap_prop = float((flat_seq.count('-')+flat_seq.count('?'))/len(flat_seq))
    n_prop = float(flat_seq.count('N')/len(flat_seq))
      
    if gap_prop >= gap_threshold:
        filter_seq = True
    elif n_prop >= n_threshold:
        filter_seq = True
    elif n_prop+gap_prop >= cumulative_threshold:
        filter_seq = True
    else:
        filter_seq = False
    
    if filter_seq:
        if write_filtered:
            file = os.path.join(os.getcwd(),"filtered_IDs.txt")
            with open(file, "a") as output_handle:
                output_handle.write(seq_name+"\n")
        else:
            return
    else:
        with open(outfile, "a") as output_handle:
            SeqIO.write(record,output_handle,"fasta")
        return
#%%
def printc(thing, level):
    '''
    Print in colour :)
    '''
    cols = {'green':'\033[1;32m', 'blue':'\033[96m'}
    col = cols[level]
    print(f"{col}{thing}\033[0m")
    return()


#%%
def main(sysargs = sys.argv[1:]):
    epilog = '''
    Analysis quits without running if outfile already exists and --force not specified.\n
    '''
    parser=argparse.ArgumentParser(description="\033[1;32mFilter Gaps and Crap\033[0m", formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    #Read inputs
    parser.add_argument('fasta', nargs="*", help="Fasta file containing sequences to check")
    parser.add_argument('-c','--cumulative-threshold',dest='cumulative_threshold', required=False, action="store", default=float(0.50), help='Maximum proportion of gaps + Ns combined')
    parser.add_argument('-g','--gap-threshold',dest='gap_threshold', required=False, action="store", default=float(0.50), help='Maximum allowed proportion of gaps in a sequence')
    parser.add_argument('-n','--n-threshold',dest='n_threshold', required=False, action="store", default=float(0.50), help='Maximum allowed proportion of Ns in a sequence')
    parser.add_argument('-o','--outfile', required=False, action="store", default=os.path.join(os.getcwd(),'alignment.filtered.fasta'), help='Name of the outfile to store results ')
    parser.add_argument('-t','--threads', required=False, action="store", default=mp.cpu_count(), help='Number of threads to use for multiprocessing of alignment')
    parser.add_argument('-w', '--write-filtered-ids',dest='write_filtered', required=False, action="store_true", default=False, help='Write filtered seq IDs to file: "filtered_IDs.txt"')
    parser.add_argument('--force', required=False, action="store_true", default=False, help='Force overwrite outfile')
    args=parser.parse_args()

    printc("\nFilter Gaps and Crap", 'green')

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    MSA = ''.join(args.fasta)

    if not os.path.exists(MSA):
        parser.print_help()
        print('#####\n\033[91mError\033[0m: Fasta file does not exist or was not specified\n#####\n')
        print('Please check your input and try again\n')
        sys.exit(-1)
    
    outfile = "alignment.filtered.fasta"
    if args.outfile:
        outfile = args.outfile
    
    if os.path.exists(outfile) and not args.force:
        parser.print_help()
        print('\n#####\n\033[91mError\033[0m: outfile already exists and overwriting not enabled\n#####\n')
        sys.exit(-1)
    elif os.path.exists(outfile) and args.force:
        print('\n#####\n\033[1;33mWarning\033[0m: outfile already exists - overwriting\n#####\n')
        os.remove(outfile)
    
    if args.write_filtered:
        if os.path.exists(os.path.join(os.getcwd(),"filtered_IDs.txt")):
            os.remove(os.path.join(os.getcwd(),"filtered_IDs.txt"))
    
    num_threads = int(mp.cpu_count())
    if args.threads:
        num_threads = int(args.threads)

    printc("\nFiltering alignment {} to remove gappy/crappy sequences.\n\nBe patient - might take a while for large files.\n\nEven if the progress bar pauses, don't panic.\n".format(''.join(args.fasta)), 'blue')
    with open(MSA) as fd, mp.Pool(num_threads) as pool:
        pool.starmap(
            filter_alignment,
            tqdm.tqdm([(record, args.gap_threshold, args.n_threshold, args.cumulative_threshold, args.outfile, args.write_filtered) for record in SeqIO.parse(fd, "fasta")]))    
    printc("\nResults in: {}\n".format(outfile), 'blue')

    sys.exit(0)

if __name__ == '__main__':
    main()

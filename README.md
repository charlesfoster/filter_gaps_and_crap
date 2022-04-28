# Filter Gaps and Crap
A simple python script to filter out gappy/crappy sequences from a multiple sequence alignment.

Author: Dr Charles Foster

# Usage
```
usage: filter_gaps_and_crap.py [-h] [-c CUMULATIVE_THRESHOLD] [-g GAP_THRESHOLD] [-n N_THRESHOLD]
                               [-o OUTFILE] [-t THREADS] [-w] [--force]
                               [fasta [fasta ...]]

Filter Gaps and Crap

positional arguments:
  fasta                 Fasta file containing sequences to check (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -c CUMULATIVE_THRESHOLD, --cumulative-threshold CUMULATIVE_THRESHOLD
                        Maximum proportion of gaps + Ns combined (default: 0.5)
  -g GAP_THRESHOLD, --gap-threshold GAP_THRESHOLD
                        Maximum allowed proportion of gaps in a sequence (default: 0.5)
  -n N_THRESHOLD, --n-threshold N_THRESHOLD
                        Maximum allowed proportion of Ns in a sequence (default: 0.5)
  -o OUTFILE, --outfile OUTFILE
                        Name of the outfile to store results
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing of alignment (default: 20)
  -w, --write-filtered-ids
                        Write filtered seq IDs to file: "filtered_IDs.txt" (default: False)
  --force               Force overwrite outfile (default: False)

Analysis quits without running if outfile already exists and --force not specified.
```

Automatically uses all cores of your CPU. Change using `--threads`.

# Dependencies
Uses default libraries, apart from the following that should be installed using (e.g.) pip:
- tqdm
- biopython

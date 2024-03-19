# pyAAscan
Python implementation of AAscan, a tool used for designing primers.
The script can be imported as a python package or used as a command line application. 
Running it from the command line requires 3 inputs:
- a text file containing the DNA sequence on which the primers should be designed (-seq)
- the index of the starting codon on the DNA sequence (-cod)
- either a single mutation (-mut) or a text file (-mutf) containing a list of mutations to make primers for

The mutation provided in -mut should be in the format [original residue][residue number][new residue], 
and the file provided in -mutf should have a single mutation per line.
Executing it from the command line will look like this:
```
python3 pyAAscan.py -seq pBAD-LEH.gb -cod 30 -mut S21R
```
Or
```
python3 pyAAscan.py -seq pBAD-LEH.gb -cod 30 -mutf mutfile.txt
```
The input file for the DNA sequence is parsed such that only lines containing numbers, spaces or ATGC characters are accepted, and only ATGC characters are saved. 
As a result, the script will work for most formats, including FASTA, EMBL, GCG and GeneBank formats. 

Various other parameters can be set, including minimum (-minl) and maximum primer length (-maxl), 
minimum (-mintm) and maximum melting temperature (-maxtm), maximum difference in melting temperature (-maxdtm),
minimum (-mino) and maximum overlap (-maxo), and minimum quality of the CGclamp (-mincg). Further info can be found using the -h or --help flag. 

AAscan's original aim as a tool for designing primers for alanine scanning can be used by invoking the --aascan flag, 
in which case providing just a residue numbers with -mut or a file containing residue numbers with -mutf is sufficient. 

Using the script as a python package allows for integration in more complex primer design workflows.
Looping over a set of desired mutations might look something like this:
```
import pyAAscan as aa

mutations = ['S21R', 'Q7P']

for mut in mutation:

    # pick the two most likely codons for a given amino acid
    codon1, codon2 = aa.BestCodons(mut[-1])

    # create the primers
    primers = aa.Mutate(seq_in=seq_in, cod1pos=30, mutpos=int(mut[1:-1]) codon1=codon1, codon2=codon2, outputmode='short')
    print(primers)

    # give the Tm of the primers independently
    for primer in primers:
        primer_noheader = primers[0].split(' ')[-1]
        print(aa.PrimerTm(primer_noheader))
```
The script requires no further packages or dependencies. 
This implementation is based on the AAscan source code provided at https://github.com/dbv123w/AAScan

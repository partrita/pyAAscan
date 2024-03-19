# pyAAscan
a python implementation of AAscan, a tool used for designing primers.
The script can either be imported as a python package or used as a command line application. 
Running it from the command line requires a minimum of 3 inputs:
- a plaintext file containing the DNA sequence on which the primers should be designed (-seq)
- the index of the starting codon on the DNA sequence (-cod)
- either a single mutation (-mut) or a file (-mutf) containing a list of mutations to make primers for
Executing it from the command line will look like this:
```
python3 pyAAscan.py -seq pBAD-LEH.gb -cod 30 -mut S21R
```
Or
```
python3 pyAAscan.py -seq pBAD-LEH.gb -cod 30 -mutf mutfile.txt
```
Various other parameters can be set, including minimum (-minl) and maximum (-maxl) primer length, 
minimum (-mintm) and maximum melting temperature (-maxtm), maximum difference in melting temperature (-maxdtm),
minimum (-mino) and maximum (-maxo) overlap, and minimum quality of the CGclamp (-mincg). further info can be found using the -h or --help flag.  
Its original aim as a tool for designing primers for alanine scanning can be used by invoking the --aascan flag, 
in which case providing just the residue numbers with -mut or -mutf is sufficient. 

Using the script as a python package allows for integration in more complex primer design workflows.
looping over a set of desired mutations might look something like this:
```
import pyAAscan as aa

mutations = ['S21R', 'Q7P']
for mut in mutation:
    # pick the two most likely codons for a given amino acid
    codon1, codon2 = aa.BestCodons(mut[-1])
    # create the primers
    primers = aa.Mutate(seq_in=seq_in, cod1pos=30, mutpos=int(mut[1:-1]) codon1=codon1, codon2=codon2)
    print(primers)
```
The script requires no further packages or dependencies. 
this implementation is based on the AAscan source code provided at https://github.com/dbv123w/AAScan

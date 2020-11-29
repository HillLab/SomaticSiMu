# SomaticSiMu
SomaticSiMu generates single and double base pair substitutions, and single base pair insertions and deletions of biologically representative mutation signature probabilities and combinations. SomaticSiMu_GUI is the GUI version of SomaticSiMu.

## Description

Simulated genomes with imposed known mutational signatures associated with cancer can be useful for benchmarking machine learning-based classifiers of genomic sequences and finetuning model hyperparameters. SomaticSiMu extracts known signature data from reference signature data, generates novel mutations on an input sequence with respect to a series of user-specified parameters, and outputs the simulated mutated sequence as a machine readable FASTA file and metadata about the position, frequency and local sequence context of each mutation. The simulation can also model temporal directed evolution across early and late stages of 37 cancer types. SomaticSiMu is developed as a lightweight, stand alone, and massively parallel software tool with a graphical user interface, built in documentation and visualization functions of mutation signature plots. The rich selection of input parameters and graphical user interface make SomaticSiMu both an easy to use application and effective as part of a wide range of experimental scenarios.  


## Installation
SomaticSiMu is implemented in Python. As long as Python is installed on your system, SomaticSiMu should run directly on your system.

$ git clone https://github.com/HillLab/SomaticSiMu\

## Base File Structure
```
├── DBS_Expected_Frequency
├── Documentation
├── Frequency_Table
├── ID_Expected_Frequency
├── Mutation_Metadata
├── Reference
├── Reference_genome
├── Sample
├── Signature_Combinations
├── kmer_ref_count
│   ├── 1-mer
│   ├── 2-mer
│   ├── 3-mer
│   ├── 4-mer
│   ├── 5-mer
│   ├── 6-mer
├── SomaticSiMu.py
├── SomaticSiMu_CC.py
├── SomaticSiMu_CC.py
```

## Quick Start

Simulate 100 sequences by imposing known mutation signatures associated with Biliary-AdenoCA onto the entire length of reference Human chromosome 22. 

```python
cd SomaticSiMu

python SomaticSiMu_GUI.py

Input Simulation Parameters: 
cancer_type = Biliary-AdenoCA
reading_frame = 1
std_outlier = 3
number_of_lineages = 100
simulation_type = end
sequence_abs_path = Homo_sapiens.GRCh38.dna.chromosome.22.fasta
slice_start = 0
slice_end = 50818467
power=1
syn_rate=1
non_syn_rate=1
```

## Parameter List

```python
"--generation", "-g", help="number of simulated sequences", default=10
"--cancer", "-c", help="cancer type"
"--reading_frame", "-f", help="index start of reading frame", default=1
"--std", "-s", help="exclude signature data outside of n std from the mean", default=3
"--simulation_type", "-v", help="simulation type", default="end"
"--slice_start", "-a", help="start of the slice of the input sequence"
"--slice_end", "-b", help="end of the slice of the input sequence"
"--power", "-p", help="multiplier of mutation burden from burden observed in in vivo samples", default=1
"--syn_rate", "-x", help="proportion of synonymous mutations out of all simulated mutations kept in the output simulated sequence", default=1
"--non_syn_rate", "-y", help="proportion of non-synonymous mutations out of all simulated mutations kept in the output simulated sequence", default=1
"--reference", "-r", help="full file path of reference sequence used as input for the simulation"
```

## Output 

Sample: Simulated sequences output into directory named after the type of cancer simulated.

Mutation_Metadata: CSV file output of each mutation simulated; the mutation type and index location on the reference input sequence. One file for each simulated sequence. 

Frequency_Table: CSV file output of summarized counts of each mutation type and local context. One file for each simulated sequence. 

Signature_Combinations: CSV file output of the signature combinations used for each iteration of the simulation. Different combinations of signatures are found operative in the same cancer type and are incorporated into the simulation. One file for each cancer type simulated. 


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)

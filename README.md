# SomaticSiMu
SomaticSiMu generates single and double base pair substitutions, and single base pair insertions and deletions of biologically representative mutation signature probabilities and combinations. SomaticSiMu_GUI is the GUI version of SomaticSiMu.

## Description

Simulated genomes with imposed known mutational signatures associated with cancer can be useful for benchmarking machine learning-based classifiers of genomic sequences and finetuning model hyperparameters. SomaticSiMu extracts known signature data from reference signature data, generates novel mutations on an input sequence with respect to a series of user-specified parameters, and outputs the simulated mutated sequence as a machine readable FASTA file and metadata about the position, frequency and local sequence context of each mutation. 

SomaticSiMu is developed as a lightweight, stand alone, and parallel software tool with a graphical user interface, built in documentation and visualization functions of mutation signature plots. The rich selection of input parameters and graphical user interface make SomaticSiMu both an easy to use application and effective as part of a wide range of experimental scenarios.  

## Requirements 

SomaticSiMu has the following dependencies:
- Python version 3.8.8 or higher
- [pandas 1.2.4 or higher](https://pypi.org/project/pandas/)
- [numpy 1.19.2 or higher](https://pypi.org/project/numpy/)
- [tqdm 4.60.0 or higher](https://pypi.org/project/tqdm/)
- [pillow 8.1.2 or higher](https://pypi.org/project/Pillow/)
- [matplotlib 3.4.1 or higher](https://pypi.org/project/matplotlib/)

Install the dependencies of SomaticSiMu to your working environment using the following command in a terminal.

```sh
pip install -r ./SomaticSiMu/requirements.txt
```

## Installation

SomaticSiMu is freely available on [GitHub](https://github.com/HillLab/SomaticSiMu). Installation requires [git](https://git-scm.com/) and [git lfs](https://git-lfs.github.com/) installed. 

Install SomaticSiMu to your working directory using the following command in a terminal.

```sh
git clone https://github.com/HillLab/SomaticSiMu\
```

## Usage

SomaticSiMu requires the absolute file path of a reference genomic sequence on the local machine as input data into the simulation. Users then select the simulation-related parameters to specify the cancer type (mutational signatures observed in whole genomes of the selected cancer type), mutation rate, location for simulated mutations, and proportion of synonymous/non-synonymous mutations as part of the simulation.

Short-Form Argument Name| Long-Form Argument Name| Argument Type | Argument Description | Argument Range 
--- | --- | --- | --- | ---
-g | --generation | Integer | Number of simulated sequences | Default = 10 ; Recommended Range: 1-100
-c | --cancer | Character | Simulated mutational signatures observed in whole genomes of the selected cancer type from PCAWG | Options: Bladder-TCC, Bone-Benign, Bone-Epith, Bone-Osteosarc, Breast-AdenoCA, Breast-DCIS, Breast-LobularCA, CNS-GBM, CNS-Medullo, CNS-Oligo, CNS-PiloAstro, Cervix-AdenoCA, Cervix-SCC, ColoRect-AdenoCA, Eso-AdenoCA, Head-SCC, Kidney-ChRCC, Kidney-RCC, Liver-HCC, Lung-AdenoCA, Lung-SCC, Lymph-BNHL, Lymph-CLL, Myeloid-AML, Myeloid-MDS, Myeloid-MPN, Ovary-AdenoCA, Panc-AdenoCA, Panc-Endocrine, Prost-AdenoCA, SKin-Melanoma, SoftTissue-Leiomyo, SoftTissue-Liposarc, Stomach-AdenoCA, Thy-AdenoCA, Uterus-AdenoCA
-f | --reading_frame | Integer | Index (1-start) of first base of the first codon in reading frame | Default = 1; Options: 1, 2, 3
-s | --std | Integer | Exclude mutational signature data from hypermutated tumors with a mutational burden `s` standard deviations from the mean mutational burden of the selected cancer type | Default = 3; Recommended Range: 0-3
-a | --slice_start | Character/Integer | Simulate mutations starting from this base index (1-start) in the reference sequence | Default = all (simulate mutations anywhere in the reference sequence), Options: Any integer from 1 up to the length of the input reference sequence
-b | --slice_end | Character/Integer | Simulate mutations starting from the slice_start index in the reference sequence up to and including this base index (1-start) |  Default = all (simulate mutations anywhere in the reference sequence), Options: Any integer greater than slice_start and up to the length of the input reference sequence
-p | --power | Integer | Multiply simulation mutation rate (baseline based on PCAWG whole genomes) by a scalar factor | Default = 1 (biologically representative) ; Recommended Range: 0.1-10
-x | --syn_rate | Float | Proportion of synonymous mutations out of all simulated mutations kept in the output simulated sequence | Default = 1 (keep all syn. mutations) ; Recommended Range: 0 (0% of syn mutations)-1 (100% of syn mutations)
-y | --non_syn_rate | Float | Proportion of non-synonymous mutations out of all simulated mutations kept in the output simulated sequence | Default = 1 (keep all non-syn. mutations) ; Recommended Range: 0 (0% of non-syn. mutations)-1 (100% of non-syn. mutations)
-r | --reference | Character | Absolute file path of reference sequence used as input for the simulation | 
-n | --normalization | Normalize mutation rates to simulate mutation types and proportions similar to the Homo Sapiens GChr38 whole genome. Different input reference sequences have different k-mer compositions compared to the whole genome that may impact the simulation of specific mutation types and their proportions.

```python
"--generation", "-g", help="number of simulated sequences", default=10
"--cancer", "-c", help="cancer type"
"--reading_frame", "-f", help="index start of reading frame", default=1
"--std", "-s", help="exclude signature data outside of n std from the mean", default=3
"--simulation_type", "-v", help="simulation type", default="end"
"--slice_start", "-a", help="start of the slice of the input sequence, default=None (start at first base)"
"--slice_end", "-b", help="end of the slice of the input sequence, default=None (end at first base)"
"--power", "-p", help="multiplier of mutation burden from burden observed in in vivo samples", default=1
"--syn_rate", "-x", help="proportion of synonymous mutations out of all simulated mutations kept in the output simulated sequence", default=1
"--non_syn_rate", "-y", help="proportion of non-synonymous mutations out of all simulated mutations kept in the output simulated sequence", default=1
"--reference", "-r", help="full file path of reference sequence used as input for the simulation"
```

## Quick Start

Difference between GUI and non-GUI version


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

## Output 

Sample: Simulated sequences output into directory named after the type of cancer simulated.

Mutation_Metadata: CSV file output of each mutation simulated; the mutation type and index location on the reference input sequence. One file for each simulated sequence. 

Frequency_Table: CSV file output of summarized counts of each mutation type and local context. One file for each simulated sequence. 

Signature_Combinations: CSV file output of the signature combinations used for each iteration of the simulation. Different combinations of signatures are found operative in the same cancer type and are incorporated into the simulation. One file for each cancer type simulated. 


## File Structure
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
```

## File Structure
Provided datasets link 

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## References

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

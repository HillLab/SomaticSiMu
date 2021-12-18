# SomaticSiMu
SomaticSiMu generates single and double base pair substitutions, and single base pair insertions and deletions of biologically representative mutation signature probabilities and combinations. SomaticSiMu_GUI is the GUI version of SomaticSiMu.

## Description

Simulated genomes with imposed known mutational signatures associated with cancer can be useful for benchmarking machine learning-based classifiers of genomic sequences and mutational signature extraction tools from mutational catalogs. SomaticSiMu extracts known signature data from reference signature data of 2,780 whole cancer genomes from the [Pan-Cancer Analysis of Whole Genomes (PCAWG) 2020 database](https://www.nature.com/articles/s41586-020-1943-3), generates novel mutations on an input reference sequence that faithfully simulate real mutational signatures, and outputs the simulated mutated DNA sequence as a machine readable FASTA file and metadata in CSV files about the position, frequency and local trinucleotide sequence context of each mutation. 

SomaticSiMu is developed as a lightweight, stand alone, and parallel software tool with an optional graphical user interface, built in documentation, and visualization functions of mutation signature plots. The rich selection of input parameters and graphical user interface make SomaticSiMu both an easy to use application and effective as part of a wide range of experimental scenarios.  

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

SomaticSiMu requires the absolute file path of a reference genomic sequence on the local machine as input data into the simulation. Users then select the simulation-related parameters (shown below) to specify the cancer type (mutational signatures observed in whole genomes of the selected cancer type), mutation rate, location for simulated mutations, and proportion of synonymous/non-synonymous mutations as part of the simulation.

The same set of arguments are offered for both `SomaticSiMu.py` and `SomaticSiMu_GUI.py`. The main difference is that `SomaticSiMu.py` is run using a terminal interface while `SomaticSiMu_GUI.py` uses a Tkinter graphical user interface to improve user accessibility along with a suite of built-in visualization functions for the simulated output data. Simulation speed and memory performance is comparable between `SomaticSiMu.py` and `SomaticSiMu_GUI.py`.

### Simulations using SomaticSiMu.py (terminal interface)

To simulate 100 unique genomic sequences using `NC_000022.11` as the reference input sequence and Skin-Melanoma associated mutational signatures, an example command in the terminal using `SomaticSiMu.py` would look like this:

```sh
python SomaticSiMu.py -g 100 -c Skin-Melanoma -r ./SomaticSiMu/Reference_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fasta
```
### Simulations using SomaticSiMu_GUI.py (graphical user interface)

To conduct the same simulation using `SomaticSiMu_GUI.py`, first run: 
```sh
python SomaticSiMu_GUI.py 
```
Then, select from drop down menus or type in the simulation parameters. Click on each simulation parameter name on the GUI interface to open up a new tab with a visual representation and description of each parameter.

### Simulation Parameters

Short-Form Argument Name| Long-Form Argument Name| Argument Type | Argument Description | Argument Range 
--- | --- | --- | --- | ---
-g | --generation | Integer | Number of simulated sequences | Default = 10 ; Recommended Range: 1-100
-c | --cancer | Character | Simulated mutational signatures observed in whole genomes of the selected cancer type from PCAWG | Options: Bladder-TCC, Bone-Benign, Bone-Epith, Bone-Osteosarc, Breast-AdenoCA, Breast-DCIS, Breast-LobularCA, CNS-GBM, CNS-Medullo, CNS-Oligo, CNS-PiloAstro, Cervix-AdenoCA, Cervix-SCC, ColoRect-AdenoCA, Eso-AdenoCA, Head-SCC, Kidney-ChRCC, Kidney-RCC, Liver-HCC, Lung-AdenoCA, Lung-SCC, Lymph-BNHL, Lymph-CLL, Myeloid-AML, Myeloid-MDS, Myeloid-MPN, Ovary-AdenoCA, Panc-AdenoCA, Panc-Endocrine, Prost-AdenoCA, SKin-Melanoma, SoftTissue-Leiomyo, SoftTissue-Liposarc, Stomach-AdenoCA, Thy-AdenoCA, Uterus-AdenoCA
-f | --reading_frame | Integer | Index (1-start) of first base of the first codon in reading frame | Default = 1; Options: 1, 2, 3
-s | --std | Integer | Exclude mutational signature data from hypermutated tumors with a mutational burden `s` standard deviations from the mean mutational burden of the selected cancer type in the PCAWG dataset | Default = 3; Recommended Range: 0-3
-a | --slice_start | Character/Integer | Simulate mutations starting from this base index (1-start) in the reference sequence | Default = all (simulate mutations anywhere in the reference sequence), Options: Any integer from 1 up to the length of the input reference sequence
-b | --slice_end | Character/Integer | Simulate mutations starting from the slice_start index in the reference sequence up to and including this base index (1-start) |  Default = all (simulate mutations anywhere in the reference sequence), Options: Any integer greater than slice_start and up to the length of the input reference sequence
-p | --power | Integer | Multiply simulation mutation rate (baseline based on PCAWG dataset) by a scalar factor | Default = 1 (biologically representative) ; Recommended Range: 0.1-10
-x | --syn_rate | Float | Proportion of synonymous mutations out of all simulated mutations kept in the output simulated sequence | Default = 1 (keep all syn. mutations) ; Recommended Range: 0 (0% of syn mutations)-1 (100% of syn mutations)
-y | --non_syn_rate | Float | Proportion of non-synonymous mutations out of all simulated mutations kept in the output simulated sequence | Default = 1 (keep all non-syn. mutations) ; Recommended Range: 0 (0% of non-syn. mutations)-1 (100% of non-syn. mutations)
-r | --reference | Character | Absolute file path of reference sequence used as input for the simulation | 
-n | --normalization | Character | Normalize mutation rates to simulate mutation types and proportions similar to the Homo Sapiens GChr38 whole genome. Different input reference sequences have different k-mer compositions compared to the whole genome that may impact the simulation of specific mutation types and their proportions. | Default: False ; Options: True, False

## Quick Start

The following quick-start examples use `SomaticSiMu.py` in a terminal interface to conduct simulation of mutational signatures. Arguments that are kept as their default value as listed in the Simulaton Parameters table are not shown for readability.

#### Simulate 10 sequences with imposed mutational signatures associated with Biliary Adenocarcinoma. Exclude hypermutants with a mutational burden that is one standard deviation beyond the mean mutational burden of the selected cancer type.

```sh
python SomaticSiMu.py -g 10 -c Biliary-AdenoCA -r ./SomaticSiMu/Reference_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fasta -s 1
```

#### Simulate 50 sequences with imposed mutational signatures associated with Colorectal Adenocarcinoma. Only simulate mutations from base index 10,000,000 to 30,000,000. 

```sh
python SomaticSiMu.py -g 50 -c ColoRect-AdenoCA -r ./SomaticSiMu/Reference_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fasta -a 10000000 -b 30000000
```

#### Simulate 100 sequences with imposed mutational signatures associated with Skin Melanoma. Increase mutation rate three-fold higher than what is observed in real tumors.

```sh
python SomaticSiMu.py -g 100 -c Skin-Melanoma -r ./SomaticSiMu/Reference_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fasta -p 3 
```

#### Simulate 30 sequences with imposed mutational signatures associated with CNS Medulloblastoma. Treat the entire input sequence as an exon. Reading frame starts at the second base of the sequence (1-start). Keep 100% of the simulated synonymous mutations. Keep 50% of the simulated non-synonymous mutations, with the other 50% randomly excluded and not present in the final output sequence.

```sh
python SomaticSiMu.py -g 30 -c CNS-Medullo -r ./SomaticSiMu/Reference_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fasta -f 2 -x 1 -y 0.5
```

#### Simulate 20 sequences with imposed mutational signatures associated with Kidney Renal Cell Carcinoma. Normalize the simulated mutations such that their types and proportions are comparable to what would be observed at the GChr38 whole genome level. 

```sh
python SomaticSiMu.py -g 20 -c Kidney-RCC -r ./SomaticSiMu/Reference_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fasta -n True
```

## Output 

Sample: Simulated sequences output into directory named after the type of cancer simulated.

Mutation_Metadata: CSV file output of each mutation simulated; the mutation type and index location on the reference input sequence. One file for each simulated sequence. 

Frequency_Table: CSV file output of summarized counts of each mutation type and local context. One file for each simulated sequence. 

Signature_Combinations: CSV file output of the signature combinations used for each iteration of the simulation. Different combinations of signatures are found operative in the same cancer type and are incorporated into the simulation. One file for each cancer type simulated. 


## File Structure
<pre>
├── Documentation                                          // Images used for documentation of SomaticSiMu-GUI.
├── Reference                                              // Reference mutational signature datasets used as the baseline for the types and proportions of mutations simulated for each cancer type.
│   └── PCAWG_sigProfiler_DBS_signatures_in_samples.csv
│   └── PCAWG_sigProfiler_ID_signatures_in_samples.csv
│   └── PCAWG_sigProfiler_SBS_signatures_in_samples.csv
│   └── sigProfiler_DBS_signatures.csv
│   └── sigProfiler_ID_signatures.csv
│   └── sigProfiler_SBS_signatures.csv
│   └── WGS_PCAWG.96.csv
│   └── WGS_PCAWG.96.xlsx
├── Reference_genome                                      // Example reference genomic sequence used for quick start simulations and testing of SomaticSiMu.
│   └── Homo_sapiens.GRCh38.dna.chromosome.22.fasta
├── kmer_ref_count                                        // 1-6mer count of GChr38 whole genome for optional normalization step during simulations.
│   └── 1-mer
│   └── 2-mer
│   └── 3-mer
│   └── 4-mer
│   └── 5-mer
│   └── 6-mer
├── Sample                                                // Output directory for simulated sequences with imposed mutations.
├── Signature_Combinations                                // Output directory that lists the set of simulated mutational signatures for each simulated sequence.
├── Mutation_Metadata                                     // Output directory that lists all simulated mutations by the mutation type, index (location in the sequence), and its local sequence context for each simulated sequence.
├── Frequency_Table                                       // Output directory that summarizes the total count of each mutation type and the possible local sequence contexts based on the SBS-96, DBS-78, or single base indel mutation classification schemes.
├── SomaticSiMu.py                                        // Simulation script (use in terminal)
├── SomaticSiMu_GUI.py                                    // Simulation script (graphic user interface)
├── requirements.txt                                      // SomaticSiMu dependencies
├── setup.py                                              // Optional automated setup script to download dependencies.
├── README.md                                             // README file.
├── LICENSE                                               // Copy of Creative Commons Attribution 4.0 International License.
</pre>

## Example Simulated Datasets
Provided datasets link 

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Citing SomaticSiMu


## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

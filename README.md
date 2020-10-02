# SomaticSiMu
SomaticSiMu generates single and double base pair substitutions, and single base pair insertions and deletions of biologically representative mutation signature probabilities and combinations. SomaticSiMu_GUI is the GUI version of SomaticSiMu.

## Installation

$ git clone https://github.com/HillLab/SomaticSiMu

## Usage

```python
cd /Users/davidchen/Documents/GitHub/SomaticSiMu

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

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)

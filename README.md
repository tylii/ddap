# DDAP
DDAP is a pathway and docking domain affinity prediction tool for type I modular PKS. 

## Compatibility
DDAP has been tested with Python 3.6.8 in Mac OS 10.14.3 and Ubuntu 18.04.1.

## Installation
DDAP's dependencies include the following tools. They can be  either installed manually or through [Bioconda](https://www.anaconda.com/distribution/).

- **(NCBI) BLAST**: [install manually](https://www.ncbi.nlm.nih.gov/books/NBK52638/), or type `conda install -c bioconda blast`
- **ClustalO**: [install manually](http://www.clustal.org/omega/), or type `conda install -c bioconda clustalo`
- **RDKit**: [install manually](https://www.rdkit.org/docs/Install.html), or type `conda install -c conda-forge rdkit`
- For Linux users: make sure you have `libxrender1`. You can install it with `sudo apt-get install libxrender1`

After downloading and unpacking the .zip file, DDAP can be installed with the following commands:

```bash
cd /path_to_ddap/ddap		# enter the root direcotry of DDAP
pip install .
``` 


## Test Run
**You can find sample input files in `/data/sample_input`**. You may run the following commands to test DDAP on sample datasets. Please rename `output_folder `. 

```bash
ddap -h						# see the help messages

# ----- Test DD Affinity Prediction ------
ddap --head ./data/sample_input/example_dd_head.fa --tail ./data/sample_input/example_dd_tail.fa -o output_folder

# ------ Test Pathway Prediction --------
ddap -i ./data/sample_input/example_pathway_input.fa -o output_folder
# alternatively
ddap -i ./data/sample_input/example_pathway_input.csv -o output_folder
ddap -i ./data/sample_input/example_KP997155.1.cluster001.gbk -o output_folder
```

## Usage

DDAP has two main functions:
**A**. predict the biosynthetic pathway of the prodcut of a type I PKS. 
**B**. predict docking domain affinity

### A. Pathway Prediction
```
ddap -i input_file -o ddap_output_folder
```
For pathway prediction, DDAP accepts `.csv`, `.fa`/`.fasta`, and antiSMASH `.gbk` files.

- `.fa`/`.fasta`: conventional sequence files of PKS peptides.

- `.csv`: the first column is the name of PKS peptides; the second culumn is the sequences of the PKS peptides; an optional third column is the monomer associated with each PKS peptide.
- `.gbk`: antiSMASH output, you can click `Download region GenBank file` in the webpage of a predicted PKS region to download the `.gbk` file. The name of the file will be similar to `NC_003888.3.region019.gbk`


### B. Docking Domain Affinity Prediction
```bash
ddap -head head.fasta -tail tail.fasta -o ddap_output_folder
``` 
Here `head.fasta` stores the C-terminal DD sequences; `tail.fasta` stores the N-terminal DD sequences. These two files must contain the same number of sequences. DDAP will predict the affinity between the frist Head and the first Tail, the second Head and the second Tail, etc.

Interpreting the Output
---------
### A. Pathway Prediction

```
ddap_output_folder/
├── Pathway_Prediction_Result/
│	  ├── pathway_prediction.csv                # Pathways Prediction Result
│	  └── structure_images/
│	  		├── pathway_1.png
│	  		├── ...
│	  		└── pathway_10.png
│	  		
└── Docking_Domain_Affinity_Prediction_Result/
    ├── DDAP_DD_prediction_result.csv     # DD Affinity Prediction Result
    ├── Model1_DD_prediction_result.csv
    ├── Model2_DD_prediction_result.csv
    │
    ├── test_files_for_BLAST/
    │   └──...
    └── test_files_for_SD_residues/
         └──...
```

- `pathway_prediction.csv` has a list of possible pathways and the associated product SMILES string. The pathways are sorted from most likely pathway and least likely pathway. The following table gives an example of the information in the `.csv` file.

| PathwayID | Pathway        | Score |
|-----------|----------------|-------|
| 1         | pksA;pksB;pksC | 0.8   |
| 2         | pksB;pksA;pksC | 0.6   |
| 3         | pksC;pksB;pksA | 0.3   |


- `structure_images/` has the product structure image for top 10 pathways.


### B. Docking Domain Affinity Prediction
```
ddap_output_folder/	  		
└── Docking_Domain_Affinity_Prediction_Result/
    ├── DDAP_DD_prediction_result.csv     # DD Affinity Prediction Result
    ├── Model1_DD_prediction_result.csv
    ├── Model2_DD_prediction_result.csv
    │
    ├── test_files_for_BLAST/
    │   └──...
    └── test_files_for_SD_residues/
         └──...
```

- `DDAP_DD_prediction_result.txt` has the predicted affinity scores of the docking domains. The following table gives an example of the information in this `.txt` file.

| Head(C-terminal-DD) | Tail(N-terminal-DD) | PredictedScore | PredictedLabel |
|---------------------|---------------------|----------------|----------------|
| head1               | tail1               | 0.4            | 0              |
| head2               | tail2               | 0.8            | 1              |
| head3               | tail3               | 0.2            | 0              |

### Other Files
The following files contains more detailed information about the prediction, including the predictors used for predictions. 

- `Model1_DD_prediction_result.csv` has the predicted scores made by Model1
- `Model2_DD_prediction_result.csv` has the predicted scores made by Model2
- `test_files_for_BLAST/` and `test_files_for_SD_residues/` have the intermediate files for Model1 and Model2.
- `DDAP_DD_prediction_result.csv `: For DD affinity prediction, this file simply lists all DD pairs and the assocaited predicted affinity. For pathway prediction, this file contains a permutaion of all possibel Head-Tail pairs, and the associated affinity. The name of the Heads are in such format `pks_H`, where pks is the name of the sequence, `_H` means this is the C-terminal DD of this sequences. Similarly, `pks_H` is the N-terminal DD of this sequence.

## DDAP databse
The data used in this study is deposited in the [DDAP databse](http://tylii.github.io/ddap).

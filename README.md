# CRISPRon-BE

## Installation

### Prerequisites

The software has been tested with the following versions, but later versions of
python, biopython, and tensorflow should also work.

* biopython  : 1.79
* python     : 3.10.8
* pandas     : 2.2.2
* tensorflow : 2.10.0
* viennarna  : 2.5.1

#### Install with Conda

The easiest way to get the prerequisites is through conda. If you have conda
installed already you can skip this step, otherwise go to
https://docs.conda.io/en/latest/miniconda.html to learn how to install conda on
your system. Once conda is correctly installed. You need to install the
CRISPRon-BE requirements with

	conda create -y -c bioconda -c conda-forge -n crispronbe \\
       python=3.10 tensorflow=2.10.0 biopython=1.79 viennarna=2.5.1 pandas=2.2.2

Later versions are also expected to work. However, the program depends on
RNAfold and versions other than 2.5.1 of the ViennaRNA package will give
slightly different results.

#### CRISPRoff, CRISPRon and Models for CRISPRon and CRISPRon-BE

CRISPRon-BE needs CRISPRoff v.1.1.2 and CRISPRon v1.0 to run in addition to the
ML models, which are quite large and are therefore separate downloads from
https://rth.dk.

The easiest way to download the models and the CRISPRon and CRISPRoff software
is to run the following script

    bin/download_and_test.sh

Alternatively you can follow the individual steps outlined in the script.

Note: CRISPRon-BE uses the CRISPRon\_V0 model, which is trained using 5
validation sets, while the models in https://github.com/RTH-tools/crispron is
trained using 6 validation sets. CRISPRon\_V0 is chosen to prevent data leakage
from predicted Cas9 indel frequency.

## Testing the software

Assuming you have installed the prerequisites in a conda environment called
crispronbe, you can run the built-in software test

	conda activate crispronbe
	./bin/test.sh

Which should end with

	Test ok

Note you are encoraged to compare your results to those found in test_org

## Running the software

Assuming you have installed the prerequisites in a conda environment called
crispronbe, you can run the software on the test data

	conda activate crispronbe
	./bin/CRISPRonBE.sh ABE test/seq.fasta test/ABE_test

To run the software on your own data, first construct a FASTA file with all the
sequences you want to have tested. See test/seq.fasta for FASTA format. Just
remember that the program needs at least 30 nucleotides to fit the full target

	prefix (4nt) -- target (20nt) -- PAM (3nt, NGG) -- suffix (3nt)

And then run the program with your own fasta file and an appropriate output
directory.

**Output files**

The script CRISPRonBE.sh outputs the following files:

- 23mers.fa: the 23 nt target + PAM sequences extracted from the input FASTA file
- 30mers.fa: the 30 nt prefix + target + PAM + suffix sequences extracted from the input FASTA file
- CRISPRparams.tsv: a tab-separated table containing the free energy changes
  computed by the CRISPRoff software. These are the RNA-DNA hybridisation
  energy (both unweighted and weighted), the DNA-DNA opening energy, the RNA-RNA
  spacer self-folding energy, and the CRISPRoff score. For details, read about
  CRISPRoff at https://rth.dk/resources/crispr/
- crispron.csv: a comma-separated table containing the 30 nt prefix + target +
  PAM + suffix sequences  and the CRISPRon predicted Cas9 indel frequencies
- crispronABE/CBE\_prediction.tsv: a tab-separated table containing the 30 nt
  prefix + target + PAM + suffix sequences, potential 30 nt edited outcomes,
  and the CRISPRon-BE predicted base editor gRNA editing efficiency, outcome
  frequency.
- log.txt: log file reporting the sequences not meeting CRISPRon-BE input requirements.

### Example run

Running on the sequence by CRISPRon-ABE

	>test
	TCAGGCTTTACAGGCCTCCGCCGCCGGGTT

You will get the following output in crispronABE\_prediction.tsv

|ID  |target                         |outcome                        |pred\_eff  |pred\_freq  |
|----|-------------------------------|-------------------------------|----------|-----------|
|test|TCAGGCTTTACAGGCCTCCGCCGCCGGGTT |TCAGGCTTTGCAGGCCTCCGCCGCCGGGTT |     50.11|      47.20|
|test|TCAGGCTTTACAGGCCTCCGCCGCCGGGTT |TCAGGCTTTGCGGGCCTCCGCCGCCGGGTT |     50.11|       2.10|
|test|TCAGGCTTTACAGGCCTCCGCCGCCGGGTT |TCAGGCTTTACGGGCCTCCGCCGCCGGGTT |     50.11|       0.81|


Here CRISPRon-ABE predicts the gRNA editing efficiency and all the potential edited outcomes.
The input is 30nt target DNA sequence, including 4nt upstream + 20nt protospacer + 3nt PAM + 3nt downstream

CRISPRon-ABE predicts this gRNA editing efficiency 49.87.
In addition, there are two \"A\"s in the editing window (from position 3 to position 10 in 20nt gRNA sequence).
After edited by ABE, this target DNA sequence may generate three different edited outcomes.
Their frequencies are 47.60, 1.92, 0.35, respectively.

### Tailormade prediction
When do prediction on new gRNAs, users can do tailormade predictions by
providing weights among training sets based on the new gRNAs most like which
datasets for the training.

The weights can be several numbers splitting by \"-\" which represents the weights for each dataset.

For ABE, 5 numbers are required, which represents the weights for SURRO-seq, Song, Arbab, Kissling ABEmax and Kissling ABE8e dataset, respectively.

For CBE, 3 numbers are required, which represents the weights for SURRO-seq, Song and Arbab dataset, respectively.

The sum of the assigned weights should be equal to 1.

We suggest the following weights for each base editor.

- ABE7.10: equal weights to both SURRO-seq and Song dataset (\"0.5-0.5-0-0-0\") since these two datasets applied a standard scaffold sequence in their experiments
- ABE8e: full weight to Kissling ABE8e (\"0-0-0-0-1\"), the only ABE8e data for training the model
- BE4: equal weights to both SURRO-seq and Song dataset (\"0.5-0.5-0\") since these two datasets applied a standard scaffold sequence in their experiments

Example of prediction using default weights

Example 1: predicting gRNAs efficiencies by ABE7.10 by default weight

    ./bin/CRISPRonBE.sh ABE test/seq.fasta test/ABE_test

    or

    ./bin/CRISPRonBE.sh ABE test/seq.fasta test/ABE_test "0.5-0.5-0-0-0"

Example of prediction using custom weights:

Example 1: predicting gRNA efficiencies measured by SURRO-seq from ABE7.10 with
weight \"1-0-0-0-0\"

    ./bin/CRISPRonBE.sh ABE test/seq.fasta test/ABE_test "1-0-0-0-0"

Example 2: providing gRNA efficiencies by assigning equal weights to Arbab and
Kissling ABE7.10 datasets in ABE7.10 with weight \"0-0-0.5-0.5-0\"

     ./bin/CRISPRonBE.sh ABE test/seq.fasta test/ABE_test "0-0-0.5-0.5-0"

## Licensing

CRISPRonBE is released under the Business Source License (BSL). It grants you
the right to copy, modify, create derivative works, redistribute, and make
non-production use of it. For production uses please
consult the Additional Use clause in the version of the BSL license provided
with the software or contact software@rth.dk for advice. Each dated version of
the license turns into the more permissive Apache License v2.0 after four
years.

Please read the complete license before usage.

## Webserver

You may wish to visit the webserver at https://rth.dk/resources/crispr, which
alows for more more analysis out-of-the-box for small scale usage. However, see
README\_webserv.md for a discussion of the differences to the off-line version
presented here.

## Citations

If you use CRISPRon-BE in your publication please cite

**Deep learning models simultaneously trained on multiple datasets improve
base-editing activity prediction.** Sun Y, Qu K, Corsi GI, Anthon C, Pan X,
Xiang X, Jensen LJ, Lin L, Luo Y, Gorodkin J. submitted.


## Contact

In case of problems or bug reports, please contact <software@rth.dk>


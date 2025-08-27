Dataset Version: 1.0
DATE: August 27th, 2025

One raw dataset generated using SURRO-seq technology Supplementary\_Data\_1.xlsx

Four processed datasets used for training and testing in CRISPRon-BE are also included in this folder.

	ABE_train.tsv.gz: training set for CRISPRon-ABE
	ABE_test.tsv.gz: test set for CRISPRon-ABE
	CBE_train.tsv.gz: training set for CRISPRon-CBE
	CBE_test.tsv.gz: test set for CRISPRon-CBE

In each file, there are 10 columns:

	refs: 20 nt gRNA sequence
	outcomes: 20 nt outcome sequence
	newc: number of reads for each outcome
	total: sum of reads from all the outcomes of one gRNA
	editingCounts: sum of reads from edited outcomes of one gRNA
	editingeff: gRNA editing efficiency, which is sum of edited reads/total reads (%) for one gRNA
	outcomefreq: outcome frequency, which is reads/totol reads (%) for one specific outcome
	surro_target: 30 nt target DNA sequence (4nt upstream + 20nt gRNA + 3nt PAM + 3nt downstream)
	source (training only): it shows the source of dataset, including "surro-seq" (SURRO-seq data), "song" (Song data), "arbab" (Arbab data), "kissling_max" (Kissling ABEmax data), "kissling_8e" (Kissling ABE8e data)
	partition: 5 partitions are used for training, and the remaining one is used for test

Note: only edited outcomes (not identical to target sequence) are used for training

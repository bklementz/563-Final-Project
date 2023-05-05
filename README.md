# 563-Final-Project
## Description:
This is the GitHub repository for my Bot563: Phylogenetic Analysis of Molecular Data semester project. 
My project is focused on resolving the internal phylogeny of the harvestman family Assamiidae (Arachnida, Opiliones, Laniatores) using sequences from 5 loci.

## Directory Layout and Contents:

notebook-log.md - contains a full description of the chosen dataset and all analyses performed over the course of the semester 
	includes: description of taxon sampling and sequenced loci; software/programs for MSA and phylogeny construction via distance, parsimony, maximum likelihood and Bayesian methods; terminal commands used in each analysis; relevant outputs.

notebook-log-final-project.md - contains descriptions of only the analyses chosen for the final project submission

./data - contains one subdirectory
	./1_Sequence_alignments - contains initial MSA fasta files for my five loci constructed by Palmieri et al. 2023, who was kind enough to lend me his dataset.

./data_clean - contains two subdirectories
	./Sequences-for-Alignment_edited - to replicate process of MSA, alignments for COI, 16S, 18S, 28S, and H3 in ./data/1_Sequence_alignments had gaps manually deleted and saved in this subdirectory.
	./Alignments - contains two subdirectories...
		./MUSCLE - contains output MSA produced for all five loci via MUSCLE
		./Mafft - contains output MSA produced for all five loci via Mafft

./manuscript - contains only the final project manuscript (563_Final-Manuscript.docx)

./analysis - contains three subdirectories, each holds all outputs produced by different phylogeny construction methods
	./ASTRAL-Outputs - contains log file of ASTRAL outputs during construction of the multiscpecies coalescence tree
		./output-speciestree-ASTRAL.log 
	./IQTree-Outputs - contains five subdirectories, one for each locus, housing all output files from IQTree Maximum Likelihood tree construction
		./16S-aligned-mafft-outputs
		./18S-aligned-mafft-outputs
		./28S-aligned-mafft-outputs
		./COI-aligned-mafft-outputs
		./H3-aligned-mafft-outputs
	./MrBayes-Outputs - contains four subdirectories, corresponding to outputs produced during Bayesian tree construction based on COI gene tree
		./COI-checkpoints - relevant outputs produced at checkpoints defined during Bayesian analysis
		./COI-run1 - outputs from the first of two runs of the Bayesian analysis
		./COI-run2 - outputs from the second of two runs of the Bayesian analysis
		./COI-total - cumulative outputs from both runs of the Bayesian analysis

./figures - contains six subdirectories
	./ASTRAL - contains two files
		./input-genetrees-ML-cat.tre - input treefile of individual gene tree topologies produced by ML for ASTRAL analysis
		./output-speciestree-ASTRAL.tre - final species tree produced by ASTRAL
	./Distance-Phylogenies - contains pdf image files of trees produced via Neighbor-Joining for 16S, COI, and H3 loci
	./IQTree-ML-Phylogenies - contains two subdirectories
		./Concatenated - contains IQ-TREE output treefiles for both concatenated analyses included in the final project 
			Sequential Concatenation analysis: concat-alignment.fasta.treefile
			Partitioned Concatenation analysis: partition.txt.treefile
		./Gene_Trees - contains IQ-TREE output treefiles for each locus
	./manuscript_figures - contains PNG files for figures included in final manuscript submission; images were produced by visualization of treefiles in FigTree v1.4.4 and exporting to Adobe Illustrator for labeling of taxa
	./MrBayes - contains the final tree topology output from MrBayes based on COI locus and a trimmed number of taxa
	./Parsimony-Phylogenies - contains pdf image files of trees produced via Maximum Parsimony for 16S, COI, and H3 loci


	
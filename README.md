# 563-Final-Project
## Description:
This is the GitHub repository for my Bot563: Phylogenetic Analysis of Molecular Data semester project. 
My project is focused on resolving the internal phylogeny of the harvestman family Assamiidae (Arachnida, Opiliones, Laniatores) using sequences from 5 loci.

## Directory Layout and Contents:
notebook-log.md - contains descriptions of: the chosen dataset (number of taxa, loci); software/programs for MSA and phylogeny construction via distance, parsimony, maximum likelihood and Bayesian methods; terminal commands used in each analysis; relevant outputs.
./data - contains one subdirectory
	./1_Sequence_alignments - contains initial MSA fasta files for my five loci constructed by Palmieri et al. 2023, who was kind enough to lend me his dataset.

./data_clean - contains two subdirectories
	./Sequences-for-Alignment_edited - to replicate process of MSA, alignments for COI, 16S, 18S, 28S, and H3 in ./data/1_Sequence_alignments had gaps manually deleted and saved in this subdirectory.
	./Alignments - contains two subdirectories...
		./MUSCLE - contains output MSA produced for all five loci via MUSCLE
		./Mafft - contains output MSA produced for all five loci via Mafft

./binaries - contains...

To-Do: edit readme here, add description of what will be contained here

./manuscript - will contain sequential drafts of the final manuscript 

./scripts - will contain extracted commands used to produce MSAs and phylogenies via distance, parsimony, maximum likelihood, and Bayesian methods.

./analysis - contains one subdirectory
	./IQTree-Outputs - contains five subdirectories, one for each locus, housing all output files from IQTree Maximum Likelihood tree construction
		./16S-aligned-mafft-outputs
		./18S-aligned-mafft-outputs
		./28S-aligned-mafft-outputs
		./COI-aligned-mafft-outputs
		./H3-aligned-mafft-outputs

./figures - contains three subdirectories
	./Distance-Phylogenies - contains pdf image files of trees produced via Neighbor-Joining for 16S, COI, and H3 loci
	./IQTree-ML-Phylogenies - contains only the tree files produced by IQTree2 for each locus
	./Parsimony-Phylogenies - contains pdf image files of trees produced via Maximum Parsimony for 16S, COI, and H3 loci


	
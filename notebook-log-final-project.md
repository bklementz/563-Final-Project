# Description of Dataset
## The Internal Phylogeny of Assamiidae (Arachnida, Opiliones, Laniatores)
Assamiidae represents a highly diverse family of Opiliones, commonly known as daddy longlegs or harvestmen, within the largely tropical suborder Laniatores.
The assamiids are characterized by a poor understanding or resolution of their internal evolutionary relationships. The traditional morphological classification system has also routinely been called into question given a plethora of monotypic genera (e.g., Mysorea) and often limited descriptions/illustrations of described taxa. 
Recently, Palmieri et al. 2023 published the first molecular phylogeny of the family, based on a 5-loci matrix for 188 taxa run via MCMCMC in MrBayes v3.2.6.
In order to assess the validity and robustness of the recovered relationships, I will reanalyze Palmieri et al.'s dataset under multispecies coalescent methods (ASTRAL), and concatenated (both with and without partitions) ML inference via IQTree

## Included Taxa
Specimens were collected across Australasia and sampled from natural history collections, including Harvard's Museum of Comparative Zoology and the Muséum d'histoire naturelle in Geneva.
Across the 5 loci, sequences are available for 199 total taxa

Relevant Outgroups
-Sequences are available for representatives of the other three suborders of Opiliones.
	-Cyphophthalmi: a single terminal representing the genus Troglosiro
	-Dyspnoi: single terminals representing the genera Hesperonemastoma, Dendrolasma, Trogulus
	-Eupnoi: single terminals representing the genera Caddo, Protolophus, and Pantopsalis
-Lantiatores Outgroups
	Infraorder: Insidiatores
		-Family Synthetonychiidae: Synthetonychia
		-Family Paranonychidae: Zuma
		-Family Travuniidae: Trojanella, Erebomaster, Theromaster, Holoscotolemon, Peltonychia
		-Family Triaenonychidae: Fumontana, Larifuga, Rostromontia, Equitius, Triaenobunus
	Infraorder: Grassatores (contains Assamiidae)
		-Family Phalangodidae: Bunofagea, Bishopella, Scotolemon (two terminals), Haasus
		-Family Sandokanidae: Martensiellus, Palaeoncopus, Sandokan, Gnomulus, Caenoncopus
		-Family Epedanidae: two unidentified genera, Pseudoepedanus, Toccolus, Nanepedanus
			-Subfamily Dibuninae: three Dibunus terminals
		-Family Zalmoxidae: Zalmoxis, Ethobunus, 
		-Family Podoctidae: Santobius (two terminals), Dongmoa, Lomanius
		-Family Stygnommatidae: Stygnomma, unidentified terminal
		-Family Samoidae: Fijicolana, Pellobunus, unidentified terminal
		-Family Biantidae: Stenostygnus (4 terminals), unidentified terminal, Metabiantes
		-Family Escadabiidae: Baculigerus, Escadabius
		-Family Guasiniidae: Guasinia
		-Family Kumulidae: Minuella, Kimula
		-Family Fissiphallius: Three terminals of Fissiphallius
		-Family Stygnopsidae: Hoplobunus, Karos, three unidentified terminals
		-Family Agoristenidae: unidentified terminal, Trinella
		-Family Nomoclastidae: Zygopachylus
		-Family Cosmetidae: Cynortula, Metalibitia, Vonones
		-Family Gonyleptidae: Goniosoma, Megapachylus, Glysterus
		-Family Stygnidae: Stygnoplus
		-Family Cranaidae: Heterocranus, Santinezia
		-Family Tithaeidae: Tithaeus
		-Family Beloniscidae: Bupares, Beloniscus
		-Family Pyramidopidae: Jarmilana, Maiorerus, Conomma
		-Incertae sedis: Arulla

Assamiidae Sampling
The remaining terminals represent nearly 95 taxa putatively identified as assamiids from throughout their range in SE Asia and the Pacific Islands
These terminals have been given an "AsXX" identifier, as well as appended with locality information. 
Other assamiids that have been secured from museum collections were often identified to existing genera, and simply identified by these genera in the terminal (e.g. Trionyxella)
Subfamily Sampling:
-Subfamily Polycoryphinae: Paktongius, As121_Laos, As110_Laos, As132_Mala, As111_Phil, As120_Thai, As109_Thai, As125_Thai, As108_Laos, As056_Gab, As021_Gab, As082_Cam, As011_Gab, As016_Gab, As010_Gab, As050_Cam, As009_Cam, As018_Gab
-Subfamily Dampetrinae: As092_EAus, As092_EDus, As104_EAus, As089_EAus, As103_NAus, As102_NAus, As101_NAus, As085_WAus, As085_WDus, As087_WAus, As095_PNG, As094_PNG, As096_PNG, As097_PNG, As134_Viet, As130_Viet, As127_Laos, As124_Indo, As133_Indo, As081_Indo, As080_Indo, As090_EAus
-Subfamily Erecinae: Neopygoplus, As116_Thai, As117_Thai, Assamiidae_DNA104858, Chilon, As034_Cam, As084_Cam, As027_Gab, As085_Cam, Assamiidae_DNA104069, Arulla_DNA102666, As114_Thai, As105_Phil, As098_Phil, As135_Viet
-Subfamily Hypoxestinae: As045_Cam, As028_Gab, As070_Lib, As083_Cam, As014_Gab, Assamiidae_DNA104857 (Bandona)
-Subfamily Mysoreinae: Assamiidae_DNA104859
-Subfamily Trionyxellinae: Trionyxella
-Subfamily Irumuinae: As032_Gab, As039_Gab, As072_Lib, As017_Gab, As031_Gab
-Subfamily Selencinae: As059_Gab, As029_Gab, As058_Gab, As040_Gab, Paraselenca, As012_Gab, As025_Gab, As057_Gab, As061_Lib, As069_Lib
Terminal identifiers can be translated to genus and species names with supplementary table S1 of Palmieri et al. (2023)

## Included Genes/Loci/Markers
DNA extraction, amplification, and sequencing yielded partial fragments of...
-two mitochondrial protein-encoding genes: 16S rRNA and Cytochrome C Oxidase subunit I
-one nuclear protein-encoding gene: histone H3
-two nuclear ribosomal genes: 18S rRNA and 28S rRNA

# Quality Control
Dataset was provided cleaned and aligned from Prashant Sharma and Luciano Palmieri-Rocha.
From Palmieri et al. 2023: sequences were assembled and cleaned using Geneious 9.1.8.
Consensus sequences were submitted to BLAST against NCBI database to identify potential contaminations.

# Alignment
-to mirror protocol of Palmieri et al. 2023, alignment will only be discussed in context of Mafft
	-for a description of alignment methods in MUSCLE, see ./563-Final-Project/notebook-log.md

## Preparation of sequences for alignment 
-sequences for each gene provided already aligned by Palmieri; included previous translation into amino acids on AliView v1.26 to remove sequences presenting premature stop codons (pseudogenes)
	-original MSA fasta files found in ./563-Final-Project/data/1_Sequence_alignments
		-relevant files = 16S.fas, 18S.fas, COI.fas, H3.fas, 28S.fas
-to personally reproduce alignment steps, gaps (designated "-" in the fasta files) were manually deleted from each of the above fasta files
Edited sequences for 16S, 18S, 28S, COI, and H3 moved to ./563-Final-Project/data_clean/Sequences-for-Alignment_Edited subdirectory


## Multisequence Alignment - Mafft

### Installation
Mafft installed via this link: <https//mafft.cbrc.jp/alignment/software/macstandard.html> and selecting the mafft-7.490-signed.pkg
Upon download, opened the .pkg file and followed installer's instructions

### Description of Algorithm(s)
Mafft offers various multiple alignment strategies
1. FFT-NS-1, FFT-NS-2 (Progressive Methods)
FFT-NS-1: first constructs a rough distance matrix by counting the number of shared 6-tuples between every sequence pair, then builds an initial tree from those distances, and finally aligns the sequences according to the branching order in the tree
	-this is the simplest option in Mafft and very fast (extremely scalable)
FFT-NS-2: since the distance matrix is only approximated at the beginning of the algorithm, this method recomputes the guide tree from the FFT-NS-1 alignment procedure and carries out a second progressive alignment from the recomputed tree

2. FFT-NS-i, NW-NS-i (Interactive Refinement Methods)
-the accuracy of the progressive alignment is improved by the iterative refinement method; implements a simplified version of PRRN
-FFT-NS-i: alignment obtained by FFT-NS-2 is subjected to iterative refinement
	-alignment is divided into two groups and realigned
	-employs tree-dependent restricted partitioning
	-process is repeated until no better scoring alignment is obtained or algorithm reaches 1000 cycles

3. L-INS-i, E-INS-i, G-INS-i (Iterative Refinement using WSP and Consistency Scores)
-Useful in difficult cases
-use new objective function combining the WSP score and COFFEE-like score to evaluate consistency between a multiple alignment and pairwise alignment
-pairwise alignments consist of three algorithms (Needleman-Wunsch, local Smith-Waterman with affine gap costs, and local alignment with generalized affine gap costs)
E-INS-i: useful in cases of conserved motifs in long unalignable regions
	-recommended if the nature of the sequences are unclear
	-assumes that the arrangement of conserved motifs are shared by all sequences
	-Unalignable residues are left unaligned at the pairwise alignment stage, because the algorithm uses the generalized affine gap cost
	-could align multiple cDNAs and multiple genomic sequences of a gene from closely related species 
	-consumes large memory space when sequences are long
L-INS-i: useful to align a set of sequences that contain subsequences flanking one alignable domain
	-flanking sequences are at first ignored in the pairwise alignment by the Smith-Waterman algorithm
	-assumes that input sequences have only one alignable domain
G-INS-i: assumes entire region can be aligned and aligns them globally with the Needleman-Wunsch algorithm
	-a set of sequences of one domain must be extracted by truncating flanking sequences

### Limitations
-accuracy: library extensions are not performed (as opposed to T-Coffee) becuase the authors believe that iterative refinement is more efficient than library extension
-scalability: if two unrelated and long genomic DNA sequences are given, FFT-NS-2 model tries to make full-length alignment and uses a large CPU time. BLAST is a more suitable method in this case
-order of alignable blocks assumed to be conserved for all input sequences

### Running the Mafft Alignments
<cd Users/bklementz/Desktop/563-Final-Project/data_clean/Sequences-for-Alignment_Edited>
-moving into the folder containing fasta files for alignment
-can use <pwd> command to ensure located in the correct subdirectory

#### Alignment of 16S Sequences:
<mafft> - command to begin alignment
'Input File?' - 16S.fasta
'Output File?' - 16S-aligned-mafft.fasta
Output format' - 4. Fasta format, input order (just need to type the number 4 and hit enter to implement choice)
Strategy - 5. L-INS-i (accurate) 
	-why L-INS-i
		-mirrors alignment strategy of Palmieri et al. 2023
		-designated the most accurate of the options included in Mafft
		-from 'https://manpages.debian.org/jessie/mafft/mafft-linsi.1': "probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information"

Program Outputs: 
command=
"/usr/local/bin/mafft"  --localpair  --maxiterate 16 --inputorder "16S.fasta" > "16S-aligned-mafft.fasta"
Type Y or just enter to run this command.
@ 
outputhat23=16
treein = 0
compacttree = 0
stacksize: 8176 kb
generating a scoring matrix for nucleotide (dist=200) ... done
All-to-all alignment.
tbfast-pair (nuc) Version 7.511
alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
0 thread(s)

outputhat23=16
Loading 'hat3.seed' ... 
done.
Writing hat3 for iterative refinement
generating a scoring matrix for nucleotide (dist=200) ... done
Gap Penalty = -1.53, +0.00, +0.00
tbutree = 1, compacttree = 0
Constructing a UPGMA tree ... 
   90 / 98
done.

Progressive alignment ... 
STEP    62 /97 
Reallocating..done. *alloclen = 1753
STEP    97 /97 
done.
tbfast (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
1 thread(s)

minimumweight = 0.000010
autosubalignment = 0.000000
nthread = 0
randomseed = 0
blosum 62 / kimura 200
poffset = 0
niter = 16
sueff_global = 0.100000
nadd = 16
Loading 'hat3' ... done.
generating a scoring matrix for nucleotide (dist=200) ... done

   90 / 98
Segment   1/  1    1- 392
STEP 004-063-1  identical.   
Converged.

done
dvtditr (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
0 thread(s)


Strategy:
 L-INS-i (Probably most accurate, very slow)
 Iterative refinement method (<16) with LOCAL pairwise alignment information
 
 Name of Final Output File: "16S-aligned-mafft.fasta" in ./data_clean/Alignments/Mafft/
 Length of final alignment (bp): 388

#### Alignment of 18S Sequences
 Same procedure as above, except...
 -input file - 18S.fasta
 -output file - 18S-aligned-mafft.fasta

Program Outputs:
command=
"/usr/local/bin/mafft"  --localpair  --maxiterate 16 --inputorder "18S.fasta" > "18S-aligned-mafft.fasta"
Type Y or just enter to run this command.
@ 
outputhat23=16
treein = 0
compacttree = 0
stacksize: 8176 kb
generating a scoring matrix for nucleotide (dist=200) ... done
All-to-all alignment.
tbfast-pair (nuc) Version 7.511
alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
0 thread(s)

outputhat23=16
Loading 'hat3.seed' ... 
done.
Writing hat3 for iterative refinement
generating a scoring matrix for nucleotide (dist=200) ... done
Gap Penalty = -1.53, +0.00, +0.00
tbutree = 1, compacttree = 0
Constructing a UPGMA tree ... 
  180 / 191
done.

Progressive alignment ... 
STEP   190 /190 
done.
tbfast (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
1 thread(s)

minimumweight = 0.000010
autosubalignment = 0.000000
nthread = 0
randomseed = 0
blosum 62 / kimura 200
poffset = 0
niter = 16
sueff_global = 0.100000
nadd = 16
Loading 'hat3' ... done.
generating a scoring matrix for nucleotide (dist=200) ... done

  180 / 191
Segment   1/  1    1-1763
STEP 002-021-1  identical.   
Converged.

done
dvtditr (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
0 thread(s)

Name of Final Output File: "18S-aligned-mafft.fasta" in ./data_clean/Alignments/Mafft/
Length of Alignment (bp): 1762

#### Alignment of 28S Sequences
 Same procedure as above, except...
 -input file - 28S.fasta
 -output file - 28S-aligned-mafft.fasta

Program Outputs:
command=
"/usr/local/bin/mafft"  --localpair  --maxiterate 16 --inputorder "28S.fasta" > "28S-aligned-mafft.fasta"
Type Y or just enter to run this command.
@ 
outputhat23=16
treein = 0
compacttree = 0
stacksize: 8176 kb
generating a scoring matrix for nucleotide (dist=200) ... done
All-to-all alignment.
tbfast-pair (nuc) Version 7.511
alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
0 thread(s)

outputhat23=16
Loading 'hat3.seed' ... 
done.
Writing hat3 for iterative refinement
generating a scoring matrix for nucleotide (dist=200) ... done
Gap Penalty = -1.53, +0.00, +0.00
tbutree = 1, compacttree = 0
Constructing a UPGMA tree ... 
  160 / 169
done.

Progressive alignment ... 
STEP   163 /168 
Reallocating..done. *alloclen = 5861
STEP   168 /168 
done.
tbfast (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
1 thread(s)

minimumweight = 0.000010
autosubalignment = 0.000000
nthread = 0
randomseed = 0
blosum 62 / kimura 200
poffset = 0
niter = 16
sueff_global = 0.100000
nadd = 16
Loading 'hat3' ... done.
generating a scoring matrix for nucleotide (dist=200) ... done

  160 / 169
Segment   1/  1    1-2461
STEP 003-057-1  identical.   
Converged.

done
dvtditr (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
0 thread(s)


Strategy:
 L-INS-i (Probably most accurate, very slow)
 Iterative refinement method (<16) with LOCAL pairwise alignment information
 
 Name of Final Output File: "28S-aligned-mafft.fasta" in ./data_clean/Alignments/Mafft
Length of Alignment (bp): 2441 

#### Alignment of COI Sequences
 Same procedure as above, except...
 -input file - COI.fasta
 -output file - COI-aligned-mafft.fasta

Program Outputs:
command=
"/usr/local/bin/mafft"  --localpair  --maxiterate 16 --inputorder "COI.fasta" > "COI-aligned-mafft.fasta"
Type Y or just enter to run this command.
@ 
outputhat23=16
treein = 0
compacttree = 0
stacksize: 8176 kb
generating a scoring matrix for nucleotide (dist=200) ... done
All-to-all alignment.
tbfast-pair (nuc) Version 7.511
alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
0 thread(s)

outputhat23=16
Loading 'hat3.seed' ... 
done.
Writing hat3 for iterative refinement
generating a scoring matrix for nucleotide (dist=200) ... done
Gap Penalty = -1.53, +0.00, +0.00
tbutree = 1, compacttree = 0
Constructing a UPGMA tree ... 
  110 / 120
done.

Progressive alignment ... 
STEP   111 /119 
Reallocating..done. *alloclen = 2317
STEP   119 /119 
done.
tbfast (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
1 thread(s)

minimumweight = 0.000010
autosubalignment = 0.000000
nthread = 0
randomseed = 0
blosum 62 / kimura 200
poffset = 0
niter = 16
sueff_global = 0.100000
nadd = 16
Loading 'hat3' ... done.
generating a scoring matrix for nucleotide (dist=200) ... done

  110 / 120
Segment   1/  1    1- 661
STEP 003-081-0  identical.   
Converged.

done
dvtditr (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
0 thread(s)


Strategy:
 L-INS-i (Probably most accurate, very slow)
 Iterative refinement method (<16) with LOCAL pairwise alignment information
 
 Name of Final Output File: "COI-aligned-mafft.fasta" in ./data_clean/Alignments/Mafft/
Length of Alignment (bp): 657

#### Alignment of H3 Sequences
Same procedure as above, except...
 -input file - H3.fasta
 -output file - H3-aligned-mafft.fasta

Program Outputs:
command=
"/usr/local/bin/mafft"  --localpair  --maxiterate 16 --inputorder "H3.fasta" > "H3-aligned-mafft.fasta"
Type Y or just enter to run this command.
@ 
outputhat23=16
treein = 0
compacttree = 0
stacksize: 8176 kb
generating a scoring matrix for nucleotide (dist=200) ... done
All-to-all alignment.
tbfast-pair (nuc) Version 7.511
alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
0 thread(s)

outputhat23=16
Loading 'hat3.seed' ... 
done.
Writing hat3 for iterative refinement
generating a scoring matrix for nucleotide (dist=200) ... done
Gap Penalty = -1.53, +0.00, +0.00
tbutree = 1, compacttree = 0
Constructing a UPGMA tree ... 
  140 / 151
done.

Progressive alignment ... 
STEP   150 /150 
done.
tbfast (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
1 thread(s)

minimumweight = 0.000010
autosubalignment = 0.000000
nthread = 0
randomseed = 0
blosum 62 / kimura 200
poffset = 0
niter = 16
sueff_global = 0.100000
nadd = 16
Loading 'hat3' ... done.
generating a scoring matrix for nucleotide (dist=200) ... done

  140 / 151
Segment   1/  1    1- 328
STEP 002-120-1  identical.   
Converged.

done
dvtditr (nuc) Version 7.511
alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
0 thread(s)


Strategy:
 L-INS-i (Probably most accurate, very slow)
 Iterative refinement method (<16) with LOCAL pairwise alignment information

Name of Final Output File: "H3-aligned-mafft.fasta" in ./data_clean/Alignments/Mafft
Length of Alignment (bp): 327

All Aligned Fasta files found in ./563-Final-Project/data_clean/Alignments/Mafft/

# Maximum Likelihood Tree Construction - Sequential and Partitioned Concatenation Analysis

## Chosen Method: IQ-Tree

### Installation
Download stable IQ-Tree version 1.6.12 from "http://www.iqtree.org/#download" (select 64-bit multicore macOS)
Upon dowload, should now have a zipped folder called "iqtree-1.6.12-MacOSX"
	found on desktop: ./Desktop/iqtree-1.6.12-MacOSX
Extracting from the zipped folder produces an accesible subdirectory iqtree-1.6.12-MacOSX
	-within this directory's subdirectory /bin, the executable application "iqtree" is functional from anywhere on the computer as long as the command paths to it (e.g. ~/Desktop/iqtree.1.6.12-MacOSX/bin/iqtree)

### Description of Algorithm
-Based on hill-climbing NNI; used to determine locally optimal trees
-NNI swaps two subtrees across internal branches
-for a given tree, likelihood scores are approximated for each NNI-tree by optimizing inner branch and adjacent 4 branches
	-only consider NNIs that increase likelihood score of current tree
-non-conflicting NNI's recorded (don't operate on same internal/adjacent branches)
-simultaneously apply all NNIs in the list to current tree and compute ML of resulting tree 
-If likelihood of resulting tree is worse, all modifications discarded except best nni - ensures a tree with better likelihood is always found
Initial Tree Generation
-Tree search start with quickly built initial tree that is improved
-IQ-Tree begins with generating 100 parsimony trees (similar to RAxML)
-the 20 trees with the best ML scores and perform hill-climbing NNI to obtain locally optimal trees
-retain 5 top topologies w/ highest likelihood in candidate tree set for optimization
Stochastic NNI 
-locally optimal trees in candidate list are randomly perturbed to escape local optima.
-then perform hill-climbing NNI on perturbed tree to reach new local optimum.
-if the new optimum has higher likelihood than the candidate tree, replaced by new tree
-tree search stops when the best tree has not changed for 100 random perturbations

### Strengths
-IQ-Tree 2 has features that reduce computational time (ultrafast bootstrap approximations, improved quartet likelihood mapping)
-IQ-Tree samples several starting trees and randomly perturbs trees during hill-climbing to prevent getting stuck in local optima-
-uses fast NNI hill-climbing algorithm, faster than traditional NNI hill climbing algorithms
-IQ-Tree2 can implement non-time reversible models, allows inference of rooted phylogenies
-given same computing time, IQ-Tree computes trees with higher likelihood than RAxML or PhyML (87.1% of benchmark datasets in Nguyen et al. 2015)
-broader exploration of tree space, given maintenance of candidate tree set which is continuously updated with better trees

### Assumptions
-samples local optima, assuming the top optima found represents the global optima/ideal tree
-IQ-Tree1 uses the same GTR (nucleotides) or WAG (amino acids) models for all tree construction
	-assumed model can only generate unrooted trees, so IQ-Tree assumes alignments are always unrooted
-I am running IQTree2 and utilizing modelfinder which allows for model selection from 286 options based on Bayesian Information Criterion scoring

### Limitations
-takes longer to run than RAxML when CPU time is variable
-IQ-Tree 1 only runs on unrooted trees
-primarily uses fast NNI, thought to produce less accurate trees than SPR
-stochastic NNI perturbations do not continue after 100 perturbations
-even with randomly sampled initial trees and perturbations, there is no guarantee the output phylogeny represents the global optimum.
-produces variation in log-likelihoods, can produce different tree topologies

### Creating Concatenated Dataset
Download phyutility
	Link: https://code.google.com/archive/p/phyutility/downloads
	Downloaded newest version: phyutility_2_2_6.tar.gz
	After download, simply opened and followed installer's instructions
	Now located on desktop: ./Desktop/phyutility
Terminal Commands
<cd ./563-Final-Project/data_clean/Alignments/Mafft> - want to work from within the subdirectory holding my individual gene alignment files
<java -jar ~/Desktop/phyutility/phyutility.jar -concat -in *.fasta -out concat-alignment.nex>
-requires pathing to phyutility.jar
*.fasta - specifies that all files ending in .fasta will be included in concatenated file

concat-alignment.nex moved to ./563-Final-Project/data_clean/Alignments/Concatenated-Alignments/

This .nex file had the header:
#NEXUS
BEGIN DATA;
	[16S-aligned-mafft.fasta_gene1 1-388 18S-aligned-mafft.fasta_gene2 389-2150 28S-aligned-mafft.fasta_gene3 2151-4591 COI-aligned-mafft.fasta_gene4 4592-5248 H3-aligned-mafft.fasta_gene5 5249-5575 ]
	DIMENSIONS NTAX=199 NCHAR=5575;
	FORMAT DATATYPE = DNA GAP = - MISSING = ?;
	MATRIX

### Running the partitioned dataset
<cd /Desktop/563-Final-Project/data_clean/Alignments/Concatenated-Alignments>
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s concat-alignment.nex>
Received error message: "ERROR: Error found at line 7, column 2125 (file position 2426):
ERROR: State specified (y) for taxon 1, character 2114, not found in list of valid symbols"

Converting concat-alignment.nex to .fasta
-opened in alignment viewer seaview, use export option and select FASTA format
	-removed partition block at the top

Creating partition block
-in new textfile:
#NEXUS 
begin sets;
	charset 16S = concat-alignment.fasta: 1-388;
	charset 18S = concat-alignment.fasta: 389-2150;
	charset 28S = concat-alignment.fasta: 2151-4591;
	charset 1COI = concat-alignment.fasta: 4592-5248\3;
	charset 2COI = concat-alignment.fasta: 4593-5248\3;
	charset 3COI = concat-alignment.fasta: 4594-5248\3;
	charset H3 = concat-alignment.fasta: 5249-5575;
end;

-multiple partitions for COI allow for variable nucleotide substitution models for each codon position, particularly useful for protein-coding genes given synonymous vs. non-synonymous mutations

file saved as partition.txt and moved to ./data_clean/Alignments/Concatenated-Alignments

Trial 2:
<~/Desktop/iqtree.1.6.12-MacOSX/bin/iqtree -s concat-alignment.fasta -q partition.txt -nt AUTO -bb 1000 -bnni -alrt 1000 -m TEST
-nt AUTO - determines optimal number of cores to use during analysis
-bb 1000 - specifies ultrafoost bootstraping with 1000 replicates
-bnni - option to reduce risk of overestimating branch supports via UF bootstrap due to model violations; ufboot will optimize each bootstrap tree using hill-climbing nearest neighbor interchange search
-alrt 1000 - performs SH-aLRT test; SH-like approximate likelihood ratio test
-m TEST implements ModelFinder to determine the best-fit model

Program outputs:
Host:    Bens-Air.lan (SSE4.2, 8 GB RAM)
Command: /Users/bklementz/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s concat-alignment.fasta -q partition.txt -nt AUTO -bb 1000 -bnni -alrt 1000 -m TEST
Seed:    8229 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Fri Apr 21 20:54:56 2023
Kernel:  SSE2 - auto-detect threads (8 CPU cores detected)

Reading partition model file partition.txt ...
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
NOTE: No CharPartition defined, use all CharSets

Loading 7 partitions...
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)
Subset	Type	Seqs	Sites	Infor	Invar	Model	Name
1		98	388	291	77	TEST	16S
2		191	1762	269	1334	TEST	18S
3		169	2441	498	1723	TEST	28S
4		120	219	133	73	TEST	1COI
5		120	219	85	102	TEST	2COI
6		120	219	219	0	TEST	3COI
7		151	327	137	170	TEST	H3
Degree of missing data: 0.175
Info: multi-threading strategy over partitions



Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.044 seconds
NOTE: ModelFinder requires 66 MB RAM!
NUMBER OF THREADS FOR PARTITION FINDING: 7
Testing GTR+F+G on supermatrix...
GTR+F+G4 / LnL: -102359.792 / df: 404 / AIC: 205527.584 / AICc: 205590.880 / BIC: 208204.508
Selecting individual models for 7 charsets using BIC...
 No. Model        Score       Charset
   1 TIM2+F+I+G4  10362.327   2COI
   2 GTR+F+I+G4   20276.832   1COI
   3 GTR+F+ASC+G4 51948.910   3COI
   4 GTR+F+I+G4   34020.738   16S
   5 GTR+F+I+G4   23570.257   H3
   6 GTR+F+I+G4   37276.230   28S
   7 TIMe+I+G4    22120.142   18S
Full partition model BIC score: 203139.153 (LnL: -99598.525  df:457)
Partition information was printed to partition.txt.best_scheme.nex
Partition information in Raxml format was printed to partition.txt.best_scheme

All model information printed to partition.txt.model.gz
CPU time for ModelFinder: 102.916 seconds (0h:1m:42s)
Wall-clock time for ModelFinder: 29.195 seconds (0h:0m:29s)
Generating 1000 samples for ultrafast bootstrap (seed: 8229)...

NOTE: 66 MB RAM (0 GB) is required!
Measuring multi-threading efficiency up to 8 CPU cores
Increase to 10 rounds for branch lengths
9 trees examined
Threads: 1 / Time: 8.140 sec / Speedup: 1.000 / Efficiency: 100% / LogL: -136516
Threads: 2 / Time: 4.630 sec / Speedup: 1.758 / Efficiency: 88% / LogL: -136516
Threads: 3 / Time: 3.729 sec / Speedup: 2.183 / Efficiency: 73% / LogL: -136516
Threads: 4 / Time: 3.625 sec / Speedup: 2.245 / Efficiency: 56% / LogL: -136516
Threads: 5 / Time: 4.276 sec / Speedup: 1.904 / Efficiency: 38% / LogL: -136516
BEST NUMBER OF THREADS: 3

Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -105508.755
Current log-likelihood at step 1: -99420.677
Current log-likelihood at step 2: -99317.376
Current log-likelihood at step 3: -99135.473
Current log-likelihood at step 4: -99098.536
Current log-likelihood at step 5: -99094.326
Current log-likelihood at step 6: -99091.739
Current log-likelihood at step 7: -99089.956
Current log-likelihood at step 8: -99089.296
Current log-likelihood at step 9: -99088.649
Current log-likelihood at step 10: -99088.270
Current log-likelihood at step 11: -99087.860
Current log-likelihood at step 12: -99087.335
Current log-likelihood at step 13: -99087.116
Current log-likelihood at step 14: -99086.926
Current log-likelihood at step 15: -99086.724
Current log-likelihood at step 16: -99086.609
Current log-likelihood at step 17: -99086.479
Current log-likelihood at step 18: -99086.358
Current log-likelihood at step 19: -99086.193
Current log-likelihood at step 20: -99086.116
Parameters optimization took 19 rounds (6.184 sec)

Computing ML distances based on estimated model parameters... 1.280 sec
WARNING: Some pairwise ML distances are too long (saturated)
Computing BIONJ tree...
0.020 seconds
Log-likelihood of BIONJ tree: -108874.111
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 6.647 second
Computing log-likelihood of 98 initial trees ... 10.430 seconds
Current best score: -99086.116

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98923.103
Current log-likelihood at step 1: -98921.863
Current log-likelihood at step 2: -98921.151
Current log-likelihood at step 3: -98920.833
Current log-likelihood at step 4: -98920.723
Current log-likelihood at step 5: -98920.506
Current log-likelihood at step 6: -98920.223
Current log-likelihood at step 7: -98920.037
Current log-likelihood at step 8: -98919.953
Parameters optimization took 7 rounds (1.944 sec)

BETTER TREE FOUND at iteration 1: -98919.953
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98902.672
Current log-likelihood at step 1: -98901.968
Current log-likelihood at step 2: -98901.798
Current log-likelihood at step 3: -98901.696
Current log-likelihood at step 4: -98901.621
Parameters optimization took 3 rounds (1.124 sec)

BETTER TREE FOUND at iteration 2: -98901.621
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98896.939
Current log-likelihood at step 1: -98896.439
Current log-likelihood at step 2: -98896.393
Parameters optimization took 1 rounds (0.517 sec)

BETTER TREE FOUND at iteration 5: -98896.393
Iteration 10 / LogL: -98910.050 / Time: 0h:1m:10s
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98841.140
Current log-likelihood at step 1: -98840.171
Current log-likelihood at step 2: -98839.993
Current log-likelihood at step 3: -98839.911
Parameters optimization took 2 rounds (0.854 sec)

BETTER TREE FOUND at iteration 11: -98839.911
Iteration 20 / LogL: -98899.311 / Time: 0h:1m:30s
Finish initializing candidate tree set (20)
Current best tree score: -98839.911 / CPU time: 57.620
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Iteration 30 / LogL: -98901.708 / Time: 0h:1m:46s (0h:4m:56s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98835.531
Current log-likelihood at step 1: -98835.387
Current log-likelihood at step 2: -98835.333
Parameters optimization took 1 rounds (0.488 sec)

BETTER TREE FOUND at iteration 31: -98835.333
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98834.965
Current log-likelihood at step 1: -98834.708
Current log-likelihood at step 2: -98834.631
Parameters optimization took 1 rounds (0.501 sec)

BETTER TREE FOUND at iteration 35: -98834.631
Iteration 40 / LogL: -98878.307 / Time: 0h:2m:3s (0h:4m:59s left)
Iteration 50 / LogL: -98836.897 / Time: 0h:2m:18s (0h:4m:0s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98834.408
Current log-likelihood at step 1: -98834.243
Current log-likelihood at step 2: -98834.201
Parameters optimization took 1 rounds (0.492 sec)

BETTER TREE FOUND at iteration 60: -98834.201
Iteration 60 / LogL: -98834.201 / Time: 0h:2m:35s (0h:4m:23s left)
Iteration 70 / LogL: -98841.723 / Time: 0h:2m:50s (0h:3m:42s left)
Iteration 80 / LogL: -98897.806 / Time: 0h:3m:5s (0h:3m:8s left)
Iteration 90 / LogL: -98900.898 / Time: 0h:3m:20s (0h:2m:37s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98832.961
Current log-likelihood at step 1: -98832.814
Current log-likelihood at step 2: -98832.796
Parameters optimization took 1 rounds (0.520 sec)

BETTER TREE FOUND at iteration 99: -98832.796
Iteration 100 / LogL: -98851.584 / Time: 0h:3m:38s (0h:3m:38s left)
Iteration 110 / LogL: -98832.800 / Time: 0h:3m:54s (0h:3m:11s left)
Iteration 120 / LogL: -98853.553 / Time: 0h:4m:9s (0h:2m:45s left)
Iteration 130 / LogL: -98861.025 / Time: 0h:4m:25s (0h:2m:21s left)
Iteration 140 / LogL: -98834.635 / Time: 0h:4m:41s (0h:1m:59s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98832.712
Current log-likelihood at step 1: -98832.689
Parameters optimization took 0 rounds (0.245 sec)

BETTER TREE FOUND at iteration 148: -98832.689
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98832.538
Current log-likelihood at step 1: -98832.504
Parameters optimization took 0 rounds (0.243 sec)

BETTER TREE FOUND at iteration 150: -98832.504
Iteration 150 / LogL: -98832.504 / Time: 0h:4m:57s (0h:3m:19s left)
Iteration 160 / LogL: -98901.168 / Time: 0h:5m:13s (0h:2m:57s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98831.836
Current log-likelihood at step 1: -98831.789
Parameters optimization took 0 rounds (0.246 sec)

BETTER TREE FOUND at iteration 166: -98831.789
Iteration 170 / LogL: -98834.620 / Time: 0h:5m:29s (0h:3m:6s left)
Iteration 180 / LogL: -98832.516 / Time: 0h:5m:45s (0h:2m:45s left)
Iteration 190 / LogL: -98832.063 / Time: 0h:6m:0s (0h:2m:24s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98831.649
Current log-likelihood at step 1: -98831.616
Parameters optimization took 0 rounds (0.230 sec)

BETTER TREE FOUND at iteration 194: -98831.616
Iteration 200 / LogL: -98837.217 / Time: 0h:6m:15s (0h:2m:57s left)
Iteration 210 / LogL: -98908.555 / Time: 0h:6m:31s (0h:2m:37s left)
Iteration 220 / LogL: -98834.658 / Time: 0h:6m:48s (0h:2m:17s left)
Iteration 230 / LogL: -98833.236 / Time: 0h:7m:4s (0h:1m:58s left)
BETTER TREE FOUND at iteration 235: -98831.609
Iteration 240 / LogL: -98845.101 / Time: 0h:7m:20s (0h:2m:55s left)
Iteration 250 / LogL: -98832.780 / Time: 0h:7m:35s (0h:2m:35s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98831.525
Current log-likelihood at step 1: -98831.513
Parameters optimization took 0 rounds (0.251 sec)

BETTER TREE FOUND at iteration 253: -98831.513
BETTER TREE FOUND at iteration 257: -98831.508
Iteration 260 / LogL: -98832.492 / Time: 0h:7m:51s (0h:2m:56s left)
Iteration 270 / LogL: -98920.438 / Time: 0h:8m:7s (0h:2m:37s left)
Iteration 280 / LogL: -98841.231 / Time: 0h:8m:22s (0h:2m:18s left)
Iteration 290 / LogL: -98843.226 / Time: 0h:8m:37s (0h:2m:0s left)
Iteration 300 / LogL: -98896.637 / Time: 0h:8m:53s (0h:1m:41s left)
Iteration 310 / LogL: -98914.221 / Time: 0h:9m:9s (0h:1m:23s left)
Iteration 320 / LogL: -98885.804 / Time: 0h:9m:25s (0h:1m:5s left)
Iteration 330 / LogL: -98842.137 / Time: 0h:9m:41s (0h:0m:47s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98831.444
Current log-likelihood at step 1: -98831.313
Current log-likelihood at step 2: -98831.297
Parameters optimization took 1 rounds (0.499 sec)

BETTER TREE FOUND at iteration 334: -98831.297
Iteration 340 / LogL: -98875.723 / Time: 0h:9m:58s (0h:2m:45s left)
Iteration 350 / LogL: -98884.352 / Time: 0h:10m:13s (0h:2m:27s left)
Iteration 360 / LogL: -98837.921 / Time: 0h:10m:29s (0h:2m:9s left)
Estimate model parameters (epsilon = 0.100)
Initial log-likelihood: -98831.261
Current log-likelihood at step 1: -98831.233
Parameters optimization took 0 rounds (0.256 sec)

BETTER TREE FOUND at iteration 361: -98831.233
Iteration 370 / LogL: -98831.891 / Time: 0h:10m:46s (0h:2m:39s left)
Iteration 380 / LogL: -98831.686 / Time: 0h:11m:1s (0h:2m:21s left)
Iteration 390 / LogL: -98840.186 / Time: 0h:11m:17s (0h:2m:3s left)
Iteration 400 / LogL: -98851.998 / Time: 0h:11m:32s (0h:1m:45s left)
Iteration 410 / LogL: -98833.629 / Time: 0h:11m:48s (0h:1m:28s left)
Iteration 420 / LogL: -98839.644 / Time: 0h:12m:5s (0h:1m:10s left)
Iteration 430 / LogL: -98864.646 / Time: 0h:12m:20s (0h:0m:53s left)
Iteration 440 / LogL: -98838.620 / Time: 0h:12m:35s (0h:0m:36s left)
Iteration 450 / LogL: -98879.070 / Time: 0h:12m:52s (0h:0m:18s left)
Iteration 460 / LogL: -98832.033 / Time: 0h:13m:9s (0h:0m:1s left)
Refining ufboot trees with NNI 5 branches...
100 samples done
200 samples done
300 samples done
400 samples done
500 samples done
600 samples done
700 samples done
800 samples done
900 samples done
1000 samples done
Total 1000 ufboot trees refined
TREE SEARCH COMPLETED AFTER 462 ITERATIONS / Time: 0h:25m:16s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
Initial log-likelihood: -98831.233
Current log-likelihood at step 1: -98831.227
Parameters optimization took 0 rounds (0.233 sec)

BEST SCORE FOUND : -98831.227

Testing tree branches by SH-like aLRT with 1000 replicates...
Alternative NNI shows better log-likelihood -98831.188 > -98831.227
4.428 sec.
Creating bootstrap support values...
Split supports printed to NEXUS file partition.txt.splits.nex
Total tree length: 23.285

Total number of iterations: 462
CPU time used for tree search: 4324.585 sec (1h:12m:4s)
Wall-clock time used for tree search: 1483.907 sec (0h:24m:43s)
Total CPU time used: 4428.490 sec (1h:13m:48s)
Total wall-clock time used: 1521.939 sec (0h:25m:21s)
UFBoot trees printed to partition.txt.ufboot

Computing bootstrap consensus tree...
Reading input file partition.txt.splits.nex...
199 taxa and 1383 splits.
Consensus tree written to partition.txt.contree
Reading input trees file partition.txt.contree
Log-likelihood of consensus tree: -98831.615

Analysis results written to: 
  IQ-TREE report:                partition.txt.iqtree
  Maximum-likelihood tree:       partition.txt.treefile
  Likelihood distances:          partition.txt.mldist
  Best partitioning scheme:      partition.txt.best_scheme.nex

Ultrafast bootstrap approximation results written to:
  Split support values:          partition.txt.splits.nex
  Consensus tree:                partition.txt.contree
  UFBoot trees:                  partition.txt.ufboot
  Screen log file:               partition.txt.log
  
partition.txt.treefile has been moved to ./figures/IQTree-ML-Phylogenies/Concatenated
remaining output files transferred to ./analysis/IQTree-Outputs/Concatenated-Outputs/Partitioned

### Running Unpartitioned Dataset
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s concat-alignment.fasta -nt AUTO -bb 1000 -bnni -alrt 1000 -m TEST>
-simply drop reference to partition file

Program Outputs:
Host:    Bens-Air.lan (SSE4.2, 8 GB RAM)
Command: /Users/bklementz/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s concat-alignment.fasta -nt AUTO -bb 1000 -bnni -alrt 1000 -m TEST
Seed:    680963 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Fri Apr 21 21:32:59 2023
Kernel:  SSE2 - auto-detect threads (8 CPU cores detected)

Reading alignment file concat-alignment.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 199 sequences with 5575 columns, 2731 distinct patterns
1632 parsimony-informative, 464 singleton sites, 3479 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                    1.51%    passed     50.21%
   2  Dongmoa                    13.13%    passed     32.06%
   3  Tithaeus                   14.62%    passed      7.62%
   4  Martensiellus               3.59%    passed     46.81%
   5  Palaeoncopus                1.51%    passed     89.94%
   6  Sandokan                    1.47%    passed     30.70%
   7  Glysterus                   2.42%    passed     58.53%
   8  Cynortula                   2.22%    passed     80.67%
   9  Zygopachylus                1.61%    passed     98.44%
  10  Santobius_DNA104931         8.09%    failed      1.55%
  11  Santobius_DNA104930        20.79%    passed     16.67%
  12  Hoplobunus                  8.75%    passed     19.38%
  13  Scotolemon_lespesi         34.03%    passed     34.88%
  14  IC_DNA104070                8.22%    failed      2.48%
  15  IC_DNA104071               19.12%    passed     21.70%
  16  Baculigerus                 5.36%    passed     37.39%
  17  Trionyxella                 2.76%    passed     81.26%
  18  Assamiidae_DNA104857       35.44%    passed     91.01%
  19  Paktongius                  7.39%    passed     62.47%
  20  Assamiidae_DNA104859       21.38%    passed     13.15%
  21  Haasus_sp                  22.40%    failed      0.49%
  22  Haasus_judaeus             35.71%    failed      2.11%
  23  Conomma                     2.30%    passed     47.44%
  24  Bishopella                 21.83%    failed      0.13%
  25  Scotolemon                  7.66%    failed      0.24%
  26  As020_Gab                  49.09%    failed      0.05%
  27  As038_Gab                  57.31%    passed      9.06%
  28  As041_Gab                  45.65%    failed      0.20%
  29  Gnomulus                    1.45%    passed     44.03%
  30  Caenoncopus                 1.22%    passed     18.62%
  31  As122_Laos                 49.26%    failed      0.01%
  32  As123_Viet                 45.31%    failed      0.01%
  33  As131_Mala                 45.58%    failed      0.69%
  34  As098_Phil                 57.74%    passed     31.80%
  35  As105_Phil                 45.90%    passed     40.04%
  36  As080_Indo                 45.72%    passed     12.79%
  37  As081_Indo                 45.60%    passed      9.70%
  38  As133_Indo                 45.11%    passed     10.06%
  39  As085_WDus                 93.27%    failed      0.00%
  40  As087_WAus                 45.45%    passed     40.33%
  41  As089_EAus                 71.91%    passed     28.44%
  42  As092_EDus                 93.27%    failed      0.00%
  43  As104_EAus                 46.82%    passed     17.48%
  44  As101_NAus                 57.52%    passed     29.13%
  45  As102_NAus                 60.39%    failed      0.00%
  46  As103_NAus                 58.40%    passed     22.23%
  47  As094_PNG                  46.44%    failed      0.76%
  48  As096_PNG                  45.63%    failed      0.33%
  49  As097_PNG                  57.33%    passed     38.03%
  50  As108_Laos                 59.57%    failed      0.01%
  51  As114_Thai                 59.59%    failed      0.00%
  52  As124_Indo                 72.02%    passed      7.68%
  53  As130_Viet                 71.32%    passed     12.14%
  54  As109_Thai                 45.17%    passed     34.88%
  55  As111_Phil                 71.59%    passed     27.78%
  56  As120_Thai                 44.97%    passed     24.98%
  57  As110_Laos                 45.42%    passed     43.32%
  58  As121_Laos                 45.09%    passed     35.50%
  59  As125_Thai                 71.77%    passed     21.92%
  60  As132_Mala                 56.86%    passed      9.12%
  61  As127_Laos                 45.51%    passed      7.64%
  62  As116_Thai                 71.41%    passed     24.19%
  63  As117_Thai                 56.93%    passed     12.09%
  64  As059_Gab                  45.61%    passed     13.29%
  65  As010_Gab                  79.93%    failed      0.19%
  66  As050_Cam                  60.02%    failed      0.00%
  67  As011_Gab                  54.04%    failed      0.00%
  68  As021_Gab                  60.84%    failed      4.42%
  69  As082_Cam                  57.43%    failed      4.30%
  70  As056_Gab                  45.58%    failed      2.95%
  71  As027_Gab                  50.89%    failed      1.21%
  72  As034_Cam                  72.41%    passed     62.23%
  73  As084_Cam                  63.75%    failed      0.16%
  74  As085_Cam                  71.93%    passed     70.54%
  75  As014_Gab                  74.08%    failed      0.05%
  76  As029_Gab                  77.90%    passed     12.08%
  77  As040_Gab                  45.43%    failed      3.50%
  78  As058_Gab                  50.80%    failed      0.16%
  79  As072_Lib                  72.14%    passed     55.15%
  80  As012_Gab                  81.27%    passed     58.01%
  81  As025_Gab                  66.31%    failed      1.04%
  82  As017_Gab                  72.93%    passed     50.53%
  83  As031_Gab                  51.96%    failed      2.87%
  84  As135_Viet                 56.90%    failed      3.62%
  85  As045_Cam                  57.33%    passed      5.90%
  86  As070_Lib                  57.52%    failed      0.54%
  87  As083_Cam                  46.22%    passed     88.08%
  88  Synthetonychia             46.96%    failed      0.00%
  89  Dendrolasma                 8.00%    failed      0.75%
  90  Trogulus                   10.71%    failed      0.00%
  91  Protolophus                 3.26%    passed     55.69%
  92  Pantopsalis                22.71%    failed      2.16%
  93  Hesperonemastoma            3.39%    failed      0.02%
  94  As061_Lib                  71.89%    passed     19.80%
  95  As069_Lib                  58.03%    failed      0.00%
  96  As032_Gab                  61.58%    failed      0.71%
  97  As039_Gab                  57.58%    failed      0.02%
  98  Troglosiro                 14.89%    failed      0.17%
  99  Pseudoepedanus             41.36%    failed      0.00%
 100  Agoristenidae_DNA105839    21.26%    passed     14.13%
 101  Trinella                   68.59%    passed     84.14%
 102  Arulla_DNA102666           17.29%    passed     57.67%
 103  Assamiidae_DNA104069       34.22%    failed      0.17%
 104  Chilon                     66.58%    failed      3.87%
 105  Paraselenca                28.22%    failed      0.00%
 106  Assamiidae_DNA104858       68.61%    passed     80.30%
 107  Fijicolana                 16.63%    passed     39.68%
 108  Stygnomma                  38.57%    passed     58.74%
 109  Pellobunus                  8.02%    passed     21.46%
 110  Icaleptes_DNA104056         8.75%    passed     59.17%
 111  Icaleptes_DNA104053        19.57%    failed      0.03%
 112  Guasinia                   31.80%    failed      0.39%
 113  Metabiantes_DNA100703      13.81%    passed     64.78%
 114  Metabiantes_DNA100704       9.15%    passed     18.39%
 115  Stenostygnus_DNA104847     13.85%    passed     89.01%
 116  Stenostygnus_DNA104849     21.99%    passed     43.64%
 117  Stenostygnus_DNA104850     21.58%    passed     14.01%
 118  Stenostygnus_DNA104848     21.99%    passed     76.80%
 119  Ethobunus                   8.50%    passed     74.18%
 120  Fissiphallius_DNA104055    19.57%    failed      0.11%
 121  Fissiphallius_DNA104057    21.87%    failed      0.07%
 122  Escadabius                 20.41%    failed      0.24%
 123  Metabiantes_DNA100335      28.09%    passed     11.03%
 124  Biantidae_DNA105668        20.25%    passed     31.02%
 125  Samoidae                   31.86%    failed      0.12%
 126  Mitraceras                 25.45%    failed      0.16%
 127  Jarmilana                  43.62%    passed     28.56%
 128  Maiorerus                  62.96%    passed     72.48%
 129  Epedanidae_DNA104861       19.61%    failed      0.02%
 130  Epedanidae_DNA104862       19.77%    failed      0.00%
 131  Goniosoma                   8.83%    passed     68.15%
 132  Megapachylus               10.82%    passed     86.25%
 133  Heterocranaus               8.02%    passed     34.06%
 134  Santinezia                  9.56%    passed     82.76%
 135  Metalibitia                19.66%    failed      0.01%
 136  Lomanius_DNA104935          8.11%    passed     60.68%
 137  Stygnopsidae_DNA103882     20.97%    passed      5.51%
 138  Stygnopsidae_DNA104855      9.52%    passed     74.65%
 139  Stygnopsidae_DNA104856      8.05%    passed     65.11%
 140  Minuella                   25.90%    failed      0.05%
 141  Karos                       8.07%    passed     94.71%
 142  Equitius                    8.05%    passed     82.67%
 143  Larifuga                   16.39%    passed     16.05%
 144  Triaenonychidae            10.40%    passed     94.65%
 145  Triaenobunus                7.95%    passed     99.13%
 146  Rostromontia               40.38%    failed      0.32%
 147  IC_DNA102668               19.66%    failed      0.09%
 148  IC_DNA102669               20.50%    failed      0.07%
 149  Zalmoxida                  14.71%    passed     36.31%
 150  IC_DNA103729               28.61%    failed      4.09%
 151  IC_DNA103572               29.65%    failed      0.09%
 152  Bunofagea                  10.78%    passed     14.13%
 153  Vonones                    77.60%    passed     11.99%
 154  Erebomaster                46.53%    failed      0.37%
 155  Holoscotolemon             19.26%    passed     29.49%
 156  Peltonychia                 7.96%    passed      8.76%
 157  Trojanella                 31.87%    passed     69.53%
 158  Epedanidae_DNA104062       28.79%    passed     69.39%
 159  Epedanidae_DNA104068       30.73%    passed     47.34%
 160  Epedanidae_DNA104066       20.43%    passed     44.19%
 161  Op106_Toccolus             42.26%    passed     16.54%
 162  Op107_Nanepedanus_rufus    42.30%    passed     16.95%
 163  Op104_Dibuninae            42.31%    passed     14.24%
 164  Op105_Dibunus              42.37%    passed     31.78%
 165  As026_Gab                  54.49%    failed      0.13%
 166  As030_Gab                  70.03%    failed      0.00%
 167  As022_Gab                  85.58%    passed     36.67%
 168  As028_Gab                  54.89%    passed     11.86%
 169  As085_WAus                 64.00%    failed      0.18%
 170  As092_EAus                 54.96%    passed     44.76%
 171  Op103_Bupares              53.78%    passed     16.94%
 172  As057_Gab                  67.46%    failed      0.00%
 173  As134_Viet                 63.50%    failed      0.05%
 174  As090_EAus                 69.31%    failed      1.27%
 175  As095_PNG                  68.04%    failed      0.00%
 176  As009_Cam                  59.09%    failed      0.19%
 177  As016_Gab                  87.26%    passed     56.00%
 178  As018_Gab                  55.35%    failed      2.71%
 179  Op049_Beloniscus           42.78%    passed     47.58%
 180  As106_Thai                 51.84%    passed      8.22%
 181  As115_Thai                 51.78%    failed      4.47%
 182  As118_Thai                 51.73%    passed      8.97%
 183  As129_Laos                 57.99%    failed      0.35%
 184  As136_Viet                 51.73%    passed     18.96%
 185  As140_Viet                 51.95%    passed     25.59%
 186  As126_Laos                 52.13%    passed     15.05%
 187  As119_Indo                 66.44%    failed      0.17%
 188  Kimula                     52.79%    failed      0.02%
 189  Fumontana                  36.61%    failed      0.14%
 190  Zuma                       55.71%    failed      0.36%
 191  Stygnoplus                 19.62%    failed      0.07%
 192  Caddo                      63.01%    passed     69.93%
 193  Theromaster                21.38%    failed      0.00%
 194  As099_Dibunus              56.72%    passed     38.37%
 195  Icaleptes                  51.01%    failed      0.00%
 196  Fissiphallius              39.37%    passed     21.40%
 197  Lomanius_DNA104934         52.07%    failed      0.00%
 198  Stygnommatidae_DNA105636   57.69%    failed      0.00%
 199  Neopygoplus                94.13%    failed      0.74%
WARNING: 77 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      39.30%  84 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.050 seconds
NOTE: ModelFinder requires 81 MB RAM!
ModelFinder will test 88 DNA models (sample size: 5575) ...
 No. Model         -LnL         df  AIC          AICc         BIC
Measuring multi-threading efficiency up to 8 CPU cores
Increase to 10 rounds for branch lengths
27 trees examined
Threads: 1 / Time: 8.203 sec / Speedup: 1.000 / Efficiency: 100% / LogL: -187133
Threads: 2 / Time: 5.168 sec / Speedup: 1.587 / Efficiency: 79% / LogL: -187133
Threads: 3 / Time: 4.332 sec / Speedup: 1.894 / Efficiency: 63% / LogL: -187133
Threads: 4 / Time: 3.672 sec / Speedup: 2.234 / Efficiency: 56% / LogL: -187133
Threads: 5 / Time: 6.627 sec / Speedup: 1.238 / Efficiency: 25% / LogL: -187133
BEST NUMBER OF THREADS: 4

  1  JC            136006.947   395 272803.893   272864.299   275421.182
  2  JC+I          117154.794   396 235101.588   235162.311   237725.502
  3  JC+G4         105827.660   396 212447.320   212508.043   215071.235
  4  JC+I+G4       104618.889   397 210031.777   210092.819   212662.318
  5  F81+F         136177.703   398 273151.406   273212.766   275788.573
  6  F81+F+I       117372.155   399 235542.310   235603.991   238186.103
  7  F81+F+G4      106117.419   399 213032.839   213094.520   215676.632
  8  F81+F+I+G4    104926.898   400 210653.796   210715.798   213304.215
  9  K2P           134154.383   396 269100.765   269161.488   271724.680
 10  K2P+I         115220.350   397 231234.700   231295.742   233865.241
 11  K2P+G4        103657.011   397 208108.022   208169.063   210738.563
 12  K2P+I+G4      102337.549   398 205471.099   205532.460   208108.266
 13  HKY+F         134334.427   399 269466.855   269528.536   272110.648
 14  HKY+F+I       115403.543   400 231607.087   231669.089   234257.506
 15  HKY+F+G4      103930.899   400 208661.797   208723.800   211312.216
 16  HKY+F+I+G4    102626.696   401 206055.392   206117.717   208712.437
 17  TNe           133649.467   397 268092.934   268153.976   270723.475
 18  TNe+I         115032.797   398 230861.593   230922.954   233498.760
 19  TNe+G4        103450.943   398 207697.885   207759.246   210335.052
 20  TNe+I+G4      102141.449   399 205080.899   205142.580   207724.692
 21  TN+F          133670.614   400 268141.229   268203.231   270791.648
 22  TN+F+I        115115.675   401 231033.350   231095.675   233690.395
 23  TN+F+G4       103577.396   401 207956.791   208019.116   210613.837
 24  TN+F+I+G4     102285.398   402 205374.795   205437.443   208038.467
 25  K3P           133973.690   397 268741.380   268802.421   271371.921
 26  K3P+I         115024.394   398 230844.788   230906.149   233481.955
 27  K3P+G4        103478.912   398 207753.823   207815.184   210390.990
 28  K3P+I+G4      102190.135   399 205178.271   205239.952   207822.064
 29  K3Pu+F        134151.193   400 269102.387   269164.389   271752.806
 30  K3Pu+F+I      115204.479   401 231210.959   231273.283   233868.004
 31  K3Pu+F+G4     103742.628   401 208287.256   208349.580   210944.301
 32  K3Pu+F+I+G4   102467.263   402 205738.527   205801.174   208402.198
 33  TPM2+F        133141.022   400 267082.044   267144.046   269732.463
 34  TPM2+F+I      114233.155   401 229268.311   229330.635   231925.356
 35  TPM2+F+G4     103273.172   401 207348.345   207410.669   210005.390
 36  TPM2+F+I+G4   102054.161   402 204912.322   204974.969   207575.993
 37  TPM2u+F       133141.022   400 267082.045   267144.047   269732.464
 38  TPM2u+F+I     114233.156   401 229268.312   229330.637   231925.357
 39  TPM2u+F+G4    103273.173   401 207348.346   207410.670   210005.391
 40  TPM2u+F+I+G4  102054.160   402 204912.321   204974.968   207575.992
 41  TPM3+F        133917.902   400 268635.805   268697.807   271286.224
 42  TPM3+F+I      114992.599   401 230787.198   230849.522   233444.243
 43  TPM3+F+G4     103728.168   401 208258.336   208320.661   210915.381
 44  TPM3+F+I+G4   102442.721   402 205689.442   205752.089   208353.113
 45  TPM3u+F       133917.902   400 268635.804   268697.806   271286.223
 46  TPM3u+F+I     114992.599   401 230787.197   230849.521   233444.242
 47  TPM3u+F+G4    103728.168   401 208258.335   208320.660   210915.381
 48  TPM3u+F+I+G4  102442.711   402 205689.422   205752.070   208353.094
 49  TIMe          133467.903   398 267731.806   267793.167   270368.973
 50  TIMe+I        114836.154   399 230470.307   230531.988   233114.100
 51  TIMe+G4       103272.665   399 207343.329   207405.010   209987.122
 52  TIMe+I+G4     101998.301   400 204796.603   204858.605   207447.022
 53  TIM+F         133486.206   401 267774.412   267836.736   270431.457
 54  TIM+F+I       114915.385   402 230634.770   230697.418   233298.441
 55  TIM+F+G4      103386.675   402 207577.349   207639.996   210241.020
 56  TIM+F+I+G4    102130.119   403 205066.239   205129.210   207736.536
 57  TIM2e         132612.115   398 266020.230   266081.591   268657.397
 58  TIM2e+I       114024.030   399 228846.060   228907.741   231489.853
 59  TIM2e+G4      102909.745   399 206617.491   206679.172   209261.284
 60  TIM2e+I+G4    101701.045   400 204202.089   204264.092   206852.508
 61  TIM2+F        132514.771   401 265831.541   265893.865   268488.586
 62  TIM2+F+I      113973.783   402 228751.566   228814.213   231415.237
 63  TIM2+F+G4     102895.218   402 206594.436   206657.084   209258.107
 64  TIM2+F+I+G4   101699.728   403 204205.456   204268.427   206875.753
 65  TIM3e         133194.948   398 267185.897   267247.258   269823.064
 66  TIM3e+I       114591.302   399 229980.604   230042.285   232624.397
 67  TIM3e+G4      103236.408   399 207270.816   207332.497   209914.608
 68  TIM3e+I+G4    101952.758   400 204705.516   204767.518   207355.935
 69  TIM3+F        133237.852   401 267277.704   267340.029   269934.750
 70  TIM3+F+I      114693.444   402 230190.888   230253.536   232854.559
 71  TIM3+F+G4     103389.667   402 207583.333   207645.980   210247.004
 72  TIM3+F+I+G4   102119.899   403 205045.798   205108.769   207716.095
 73  TVMe          132630.071   399 266058.142   266119.824   268701.935
 74  TVMe+I        113730.169   400 228260.338   228322.340   230910.757
 75  TVMe+G4       102857.740   400 206515.480   206577.482   209165.899
 76  TVMe+I+G4     101646.010   401 204094.020   204156.345   206751.065
 77  TVM+F         132707.771   402 266219.542   266282.189   268883.213
 78  TVM+F+I       113793.988   403 228393.975   228456.946   231064.272
 79  TVM+F+G4      103000.597   403 206807.195   206870.166   209477.492
 80  TVM+F+I+G4    101802.252   404 204412.504   204475.800   207089.427
 81  SYM           132141.036   400 265082.073   265144.075   267732.492
 82  SYM+I         113555.406   401 227912.812   227975.137   230569.857
 83  SYM+G4        102630.544   401 206063.088   206125.412   208720.133
 84  SYM+I+G4      101435.917   402 203675.833   203738.481   206339.505
 85  GTR+F         132064.812   403 264935.623   264998.594   267605.920
 86  GTR+F+I       113523.272   404 227854.544   227917.840   230531.467
 87  GTR+F+G4      102638.793   404 206085.587   206148.883   208762.510
 88  GTR+F+I+G4    101455.321   405 203720.641   203784.263   206404.191
Akaike Information Criterion:           SYM+I+G4
Corrected Akaike Information Criterion: SYM+I+G4
Bayesian Information Criterion:         SYM+I+G4
Best-fit model: SYM+I+G4 chosen according to BIC

All model information printed to concat-alignment.fasta.model.gz
CPU time for ModelFinder: 294.939 seconds (0h:4m:54s)
Wall-clock time for ModelFinder: 82.523 seconds (0h:1m:22s)
Generating 1000 samples for ultrafast bootstrap (seed: 680963)...

NOTE: 81 MB RAM (0 GB) is required!
Measuring multi-threading efficiency up to 8 CPU cores
Increase to 10 rounds for branch lengths
10 trees examined
Threads: 1 / Time: 8.573 sec / Speedup: 1.000 / Efficiency: 100% / LogL: -136892
Threads: 2 / Time: 5.236 sec / Speedup: 1.637 / Efficiency: 82% / LogL: -136892
Threads: 3 / Time: 3.784 sec / Speedup: 2.265 / Efficiency: 76% / LogL: -136892
Threads: 4 / Time: 3.229 sec / Speedup: 2.655 / Efficiency: 66% / LogL: -136892
Threads: 5 / Time: 6.497 sec / Speedup: 1.319 / Efficiency: 26% / LogL: -136892
BEST NUMBER OF THREADS: 4

Estimate model parameters (epsilon = 0.100)
Thoroughly optimizing +I+G parameters from 10 start values...
Init pinv, alpha: 0.000, 1.000 / Estimate: 0.000, 0.212 / LogL: -102630.484
Init pinv, alpha: 0.069, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.919
Init pinv, alpha: 0.139, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.886
Init pinv, alpha: 0.208, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.873
Init pinv, alpha: 0.277, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.874
Init pinv, alpha: 0.347, 1.000 / Estimate: 0.588, 0.624 / LogL: -101435.868
Init pinv, alpha: 0.416, 1.000 / Estimate: 0.588, 0.624 / LogL: -101435.870
Init pinv, alpha: 0.485, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.865
Init pinv, alpha: 0.555, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.897
Init pinv, alpha: 0.624, 1.000 / Estimate: 0.588, 0.625 / LogL: -101435.952
Optimal pinv,alpha: 0.588, 0.625 / LogL: -101435.865

Parameters optimization took 27.147 sec
Computing ML distances based on estimated model parameters... 1.307 sec
WARNING: Some pairwise ML distances are too long (saturated)
Computing BIONJ tree...
0.019 seconds
Log-likelihood of BIONJ tree: -111283.974
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 5.728 second
Computing log-likelihood of 98 initial trees ... 8.078 seconds
Current best score: -101156.315

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -100983.392
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 2: -100976.704
Iteration 10 / LogL: -101023.602 / Time: 0h:1m:24s
Iteration 20 / LogL: -101003.097 / Time: 0h:1m:38s
Finish initializing candidate tree set (20)
Current best tree score: -100976.704 / CPU time: 43.635
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 21: -100973.759
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 30: -100966.195
Iteration 30 / LogL: -100966.195 / Time: 0h:1m:51s (0h:6m:23s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 34: -100965.103
Iteration 40 / LogL: -100970.606 / Time: 0h:2m:2s (0h:4m:56s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 41: -100946.390
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 48: -100943.966
Iteration 50 / LogL: -100949.268 / Time: 0h:2m:16s (0h:4m:32s left)
Iteration 60 / LogL: -100956.247 / Time: 0h:2m:30s (0h:3m:43s left)
Iteration 70 / LogL: -101027.066 / Time: 0h:2m:43s (0h:3m:4s left)
Iteration 80 / LogL: -100946.380 / Time: 0h:2m:55s (0h:2m:31s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 82: -100942.264
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 83: -100939.760
Iteration 90 / LogL: -100944.698 / Time: 0h:3m:8s (0h:3m:17s left)
Iteration 100 / LogL: -100940.179 / Time: 0h:3m:20s (0h:2m:48s left)
Iteration 110 / LogL: -100942.624 / Time: 0h:3m:32s (0h:2m:22s left)
Iteration 120 / LogL: -100949.092 / Time: 0h:3m:46s (0h:1m:59s left)
Iteration 130 / LogL: -100974.469 / Time: 0h:3m:58s (0h:1m:37s left)
Iteration 140 / LogL: -100953.219 / Time: 0h:4m:11s (0h:1m:17s left)
Iteration 150 / LogL: -100944.218 / Time: 0h:4m:23s (0h:0m:58s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 157: -100939.243
Iteration 160 / LogL: -100939.793 / Time: 0h:4m:35s (0h:2m:48s left)
Iteration 170 / LogL: -100940.298 / Time: 0h:4m:48s (0h:2m:28s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 177: -100939.148
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 179: -100938.843
Iteration 180 / LogL: -100943.600 / Time: 0h:5m:1s (0h:2m:46s left)
Iteration 190 / LogL: -100950.651 / Time: 0h:5m:13s (0h:2m:27s left)
Iteration 200 / LogL: -100956.584 / Time: 0h:5m:25s (0h:2m:9s left)
Estimate model parameters (epsilon = 0.100)
UPDATE BEST LOG-LIKELIHOOD: -100938.558
Iteration 210 / LogL: -101037.164 / Time: 0h:5m:37s (0h:1m:51s left)
Iteration 220 / LogL: -100943.545 / Time: 0h:5m:50s (0h:1m:34s left)
Iteration 230 / LogL: -100976.831 / Time: 0h:6m:3s (0h:1m:17s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 240: -100938.282
Iteration 240 / LogL: -100938.282 / Time: 0h:6m:15s (0h:2m:37s left)
Iteration 250 / LogL: -100947.343 / Time: 0h:6m:27s (0h:2m:20s left)
Estimate model parameters (epsilon = 0.100)
UPDATE BEST LOG-LIKELIHOOD: -100937.977
Iteration 260 / LogL: -100941.392 / Time: 0h:6m:40s (0h:2m:3s left)
Iteration 270 / LogL: -100939.899 / Time: 0h:6m:53s (0h:1m:47s left)
Iteration 280 / LogL: -100942.388 / Time: 0h:7m:6s (0h:1m:31s left)
Iteration 290 / LogL: -100938.675 / Time: 0h:7m:18s (0h:1m:15s left)
Iteration 300 / LogL: -100947.794 / Time: 0h:7m:30s (0h:1m:0s left)
Iteration 310 / LogL: -100938.555 / Time: 0h:7m:43s (0h:0m:44s left)
Iteration 320 / LogL: -100938.801 / Time: 0h:7m:54s (0h:0m:29s left)
Iteration 330 / LogL: -100943.265 / Time: 0h:8m:6s (0h:0m:14s left)
Iteration 340 / LogL: -100941.250 / Time: 0h:8m:18s (0h:0m:0s left)
Refining ufboot trees with NNI 1 branch...
100 samples done
200 samples done
300 samples done
400 samples done
500 samples done
600 samples done
700 samples done
800 samples done
900 samples done
1000 samples done
Total 1000 ufboot trees refined
TREE SEARCH COMPLETED AFTER 341 ITERATIONS / Time: 0h:12m:7s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -100937.977
Optimal log-likelihood: -100937.968
Rate parameters:  A-C: 1.24180  A-G: 3.35693  A-T: 2.99290  C-G: 0.81424  C-T: 5.85773  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Proportion of invariable sites: 0.589
Gamma shape alpha: 0.624
Parameters optimization took 1 rounds (0.215 sec)
BEST SCORE FOUND : -100937.968

Testing tree branches by SH-like aLRT with 1000 replicates...
3.577 sec.
Creating bootstrap support values...
Split supports printed to NEXUS file concat-alignment.fasta.splits.nex
Total tree length: 12.222

Total number of iterations: 341
CPU time used for tree search: 2628.763 sec (0h:43m:48s)
Wall-clock time used for tree search: 671.777 sec (0h:11m:11s)
Total CPU time used: 2831.208 sec (0h:47m:11s)
Total wall-clock time used: 731.872 sec (0h:12m:11s)
UFBoot trees printed to concat-alignment.fasta.ufboot

Computing bootstrap consensus tree...
Reading input file concat-alignment.fasta.splits.nex...
199 taxa and 1244 splits.
Consensus tree written to concat-alignment.fasta.contree
Reading input trees file concat-alignment.fasta.contree
Log-likelihood of consensus tree: -100938.996

Analysis results written to: 
  IQ-TREE report:                concat-alignment.fasta.iqtree
  Maximum-likelihood tree:       concat-alignment.fasta.treefile
  Likelihood distances:          concat-alignment.fasta.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          concat-alignment.fasta.splits.nex
  Consensus tree:                concat-alignment.fasta.contree
  UFBoot trees:                  concat-alignment.fasta.ufboot
  Screen log file:               concat-alignment.fasta.log

concat-alignment.fasta.treefile transferred to ./figures/IQTree-ML-Phylogenies/Concatenated
remaining output files transferred to ./analysis/IQTree-Outputs/Concatenated-Outputs/Unpartitioned

# Coalescent Methods
## Chosen Method: ASTRAL
### Installation
From "https://github.com/smirarab/ASTRAL/blob/master/README.md#installation",
	-simply select the zip file link and extract the contents
	-ASTRAL folder relocated to desktop - ./Desktop/Astral/
### Description of Algorithm
### Description of Algorithm
-A method for constructing a species tree from a set of gene trees
-functions under the multi-species coalescent model of gene evolution
	-incorporates incomplete lineage sorting and disagreement between gene trees and the species tree (differing topologies)
	-species tree is inferred from a distribution of gene tree topologies, the species tree topology is the most probable gene tree topology
		-can count the number of gene trees and pick the most frequent one as the species tree
-one of the scalable alternatives to co-estimation of gene and species trees under a joint inference
-ASTRAL is a "summary method"
	-these methods first infer gene trees independently for all loci, then combine the gene trees to form a species tree
	-effectively ignores the dependent nature of different loci
-Input: a set of unrooted gene trees, ASTRAL will find the species tree that agrees with the largest number of quartet trees from the set of gene trees
	-Heuristic version: used on large datasets, constrains search space (reduce run time)
		-includes a set of bipartitions as part of input, and output species tree must draw its bipartitions from that set
		-reduces the NP-hard problem to be solved in polynomial time (Mirarab et al. 2014)
		-also called the Maximum Quartet Support Species Tree approach (MQSST)
		-MQSST uses dynamic programming to solve optimization problem
		-Use set of defined bipartitions to generate series of tripartitions
		-for each tripartition, calculate the number of quartet trees induced by input gene trees associated to that tripartition. 
		-species tree constructed by calculating a score for individual tripartitions based on recursive formula
-A binary tree can be represented as a set of tripartitions, one per node
-DP can be used to find a set of tripartitions that can be combined into a binary tree that have the maximum number of possible shared quartets with gene trees
### Strengths
-more accurate than concatenated ML phylogenies due to higher statistical consistency
-breaking analysis into independent steps increases scalability
-MQSST accounts for estimation error of dominant quartet tree 
	-takes into account relative frequency of all three quartet topologies and weights them; other methods try to first find dominant quartet topology and summarize them
-statistically consistent when input gene trees are sampled randomly under MSC model (no sampling bias/model violations)
-remains consistent when species are missing from gene trees
-more accurate than concatenation with ML when gene tree discordance is high (high ILS)
-tends to outperform NJst by small margins in simulations (Mirarab and Warnow 2015)
-dominates MP-EST when using large datasets especially

### Limitations
-ASTRAL are statistially inconsistent if each gene has limited length and gene trees constructed with maximum likelihood
	-fails due to LBA
-statistically inconsistent if gene trees evolving in a network (ILS+gene flow)
-concatenated ML more accurate when gene tree discordance is low, or when gene tree error is high
-some scalability concerns, "increasing the number of species should roughly quadruple the running time" - (Mirarab 2019; class paper)
### Assumptions
-robustness to missing data in exact version requires that presence of a gene for a species is independent of of gene tree topology and presence of other genes for that species
	-heuristic version requires that each clade in species tree has a non-zero chance of having no missing data in each gene

### Creating Individual Gene Trees in IQ-TREE

<cd ~/Desktop/563-Final-Project/data/Alignments/Mafft>
-navigating to directory containing my MSAs

#### 16S Gene Tree
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s 16S-aligned-mafft.fasta>
-running 16S Mafft alignment 
-first specify location of iqtree application in downloaded directory
then -s <name of input file>
Output: 
Reading alignment file 16S-aligned-mafft.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
Alignment has 98 sequences with 388 columns, 335 distinct patterns
291 parsimony-informative, 20 singleton sites, 77 constant sites
                      Gap/Ambiguity  Composition  p-value
   1  Zalmoxis                6.44%    passed      6.92%
   2  Dongmoa                 6.19%    failed      0.28%
   3  Tithaeus                6.19%    failed      1.27%
   4  Martensiellus           6.44%    passed     24.08%
   5  Palaeoncopus            6.44%    passed     39.49%
   6  Sandokan                6.44%    failed      0.78%
   7  Glysterus               6.19%    failed      0.41%
   8  Cynortula               6.96%    passed     10.81%
   9  Zygopachylus            6.19%    failed      4.50%
  10  Santobius_DNA104931     6.19%    failed      4.79%
  11  Santobius_DNA104930    20.88%    passed     13.88%
  12  Hoplobunus              6.44%    passed     10.08%
  13  Scotolemon_lespesi      7.73%    failed      0.25%
  14  IC_DNA104070            7.73%    failed      2.16%
  15  IC_DNA104071            6.44%    passed     20.84%
  16  Baculigerus             7.47%    failed      0.30%
  17  Trionyxella             6.44%    passed     45.42%
  18  Assamiidae_DNA104857    6.96%    passed     96.50%
  19  Paktongius              6.19%    failed      3.98%
  20  Assamiidae_DNA104859    6.19%    passed     15.16%
  21  Haasus_sp               9.02%    passed     13.45%
  22  Haasus_judaeus         15.98%    passed     14.95%
  23  Conomma                18.04%    failed      4.66%
  24  Bishopella              6.70%    failed      0.83%
  25  Scotolemon              7.99%    passed     13.08%
  26  As020_Gab               4.12%    failed      0.18%
  27  As038_Gab               3.61%    failed      0.12%
  28  As041_Gab               3.61%    failed      0.40%
  29  Gnomulus                3.61%    failed      3.18%
  30  Caenoncopus             3.35%    failed      3.98%
  31  As122_Laos              3.61%    passed     25.72%
  32  As123_Viet              3.61%    failed      4.83%
  33  As131_Mala              3.61%    passed     42.41%
  34  As098_Phil              3.35%    passed      8.23%
  35  As105_Phil              3.35%    passed     32.06%
  36  As080_Indo              4.12%    passed     88.81%
  37  As081_Indo              3.87%    passed     51.98%
  38  As133_Indo              3.61%    passed     59.21%
  39  As085_WDus              3.35%    passed     75.49%
  40  As087_WAus              3.61%    passed     73.09%
  41  As089_EAus              3.35%    passed     66.99%
  42  As092_EDus              3.35%    passed     95.88%
  43  As104_EAus              3.35%    passed     86.25%
  44  As101_NAus              3.35%    passed     87.05%
  45  As102_NAus              3.35%    passed     80.51%
  46  As103_NAus              3.35%    passed     86.15%
  47  As094_PNG               3.61%    passed     61.01%
  48  As096_PNG               3.35%    passed     48.47%
  49  As097_PNG               3.35%    passed     72.40%
  50  As108_Laos              3.61%    passed     80.57%
  51  As114_Thai              3.35%    passed     56.13%
  52  As124_Indo              4.12%    passed     12.63%
  53  As130_Viet              3.35%    passed     10.98%
  54  As109_Thai              3.87%    passed     38.36%
  55  As111_Phil              3.87%    passed     78.54%
  56  As120_Thai              3.87%    passed     53.07%
  57  As110_Laos              3.35%    passed     43.38%
  58  As121_Laos              3.35%    passed     60.90%
  59  As125_Thai              8.51%    passed     72.55%
  60  As132_Mala              3.35%    passed     50.46%
  61  As127_Laos              4.38%    failed      0.55%
  62  As116_Thai              3.87%    passed     87.05%
  63  As117_Thai              4.12%    passed     52.01%
  64  As059_Gab               3.61%    passed     48.10%
  65  As010_Gab               3.35%    passed      5.39%
  66  As050_Cam               3.35%    failed      3.63%
  67  As011_Gab               3.35%    passed     12.22%
  68  As021_Gab               3.35%    failed      4.39%
  69  As082_Cam               3.35%    failed      3.73%
  70  As056_Gab               3.35%    passed      6.69%
  71  As027_Gab               3.35%    failed      0.59%
  72  As034_Cam               3.61%    failed      0.27%
  73  As084_Cam              48.20%    passed      7.57%
  74  As085_Cam               3.35%    failed      0.08%
  75  As014_Gab               3.35%    passed     14.03%
  76  As029_Gab               3.35%    failed      0.13%
  77  As040_Gab               3.35%    failed      0.01%
  78  As058_Gab               3.35%    failed      0.04%
  79  As072_Lib               3.87%    failed      0.01%
  80  As012_Gab               3.35%    failed      0.00%
  81  As025_Gab               4.12%    failed      0.00%
  82  As017_Gab               3.35%    failed      0.91%
  83  As031_Gab               3.61%    failed      0.71%
  84  As135_Viet              3.35%    failed      0.05%
  85  As045_Cam               3.87%    passed     65.35%
  86  As070_Lib               3.35%    failed      0.04%
  87  As083_Cam               3.35%    failed      2.03%
  88  Synthetonychia          7.73%    failed      0.17%
  89  Dendrolasma             7.47%    passed     82.74%
  90  Trogulus                7.73%    passed      7.53%
  91  Protolophus             7.22%    passed     52.38%
  92  Pantopsalis             9.28%    passed      8.08%
  93  Hesperonemastoma       11.34%    passed     21.95%
  94  As061_Lib               3.35%    failed      0.00%
  95  As069_Lib               3.61%    failed      0.00%
  96  As032_Gab               5.41%    failed      0.00%
  97  As039_Gab               5.41%    failed      0.00%
  98  Troglosiro              8.51%    failed      0.39%
****  TOTAL                   5.62%  42 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.007 seconds
NOTE: ModelFinder requires 10 MB RAM!
ModelFinder will test 286 DNA models (sample size: 388) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  JC            18869.039    193 38124.077    38510.077    38888.551
  2  JC+I          18006.162    194 36400.323    36792.344    37168.758
  3  JC+G4         17142.697    194 34673.393    35065.414    35441.828
  4  JC+I+G4       17112.067    195 34614.134    35012.259    35386.531
  5  JC+R2         17343.890    195 35077.780    35475.905    35850.176
  6  JC+R3         17174.213    197 34742.427    35153.016    35522.745
  7  JC+R4         17097.876    199 34593.752    35017.156    35381.992
  8  JC+R5         17088.060    201 34578.121    35014.701    35374.283
  9  JC+R6         17083.422    203 34572.845    35022.975    35376.929
 14  F81+F         18794.432    196 37980.863    38385.178    38757.220
 15  F81+F+I       17896.289    197 36186.579    36597.168    36966.897
 16  F81+F+G4      17000.582    197 34395.165    34805.754    35175.483
 17  F81+F+I+G4    16964.813    198 34325.626    34742.578    35109.905
 18  F81+F+R2      17192.706    198 34781.412    35198.364    35565.691
 19  F81+F+R3      17018.677    200 34437.354    34867.301    35229.556
 20  F81+F+R4      16945.659    202 34295.319    34738.627    35095.442
 21  F81+F+R5      16932.910    204 34273.821    34730.870    35081.866
 22  F81+F+R6      16927.758    206 34267.517    34738.699    35083.484
 27  K2P           18532.086    194 37452.173    37844.194    38220.608
 28  K2P+I         17660.794    195 35711.587    36109.712    36483.983
 29  K2P+G4        16765.316    195 33920.632    34318.757    34693.029
 30  K2P+I+G4      16728.954    196 33849.908    34254.222    34626.265
 31  K2P+R2        16943.636    196 34279.271    34683.586    35055.628
 32  K2P+R3        16770.388    198 33936.776    34353.729    34721.055
 33  K2P+R4        16713.933    200 33827.866    34257.813    34620.067
 34  K2P+R5        16697.156    202 33798.311    34241.620    34598.435
 35  K2P+R6        16691.185    204 33790.370    34247.420    34598.416
 36  K2P+R7        16691.068    206 33794.136    34265.318    34610.103
 40  HKY+F         18384.556    197 37163.111    37573.701    37943.429
 41  HKY+F+I       17449.579    198 35295.157    35712.110    36079.436
 42  HKY+F+G4      16492.733    198 33381.465    33798.418    34165.744
 43  HKY+F+I+G4    16441.878    199 33281.755    33705.159    34069.995
 44  HKY+F+R2      16691.967    199 33781.933    34205.338    34570.173
 45  HKY+F+R3      16505.683    201 33413.366    33849.947    34209.528
 46  HKY+F+R4      16429.929    203 33265.858    33715.989    34069.942
 47  HKY+F+R5      16412.964    205 33235.928    33699.993    34047.934
 48  HKY+F+R6      16408.508    207 33231.016    33709.416    34050.945
 53  TNe           18521.144    195 37432.287    37830.412    38204.683
 54  TNe+I         17646.709    196 35685.419    36089.733    36461.776
 55  TNe+G4        16755.088    196 33902.175    34306.490    34678.532
 56  TNe+I+G4      16720.412    197 33834.824    34245.413    34615.142
 57  TNe+R2        16922.427    197 34238.854    34649.444    35019.172
 58  TNe+R3        16755.100    199 33908.199    34331.604    34696.439
 59  TNe+R4        16703.331    201 33808.662    34245.243    34604.824
 60  TNe+R5        16685.743    203 33777.486    34227.616    34581.570
 61  TNe+R6        16681.867    205 33773.733    34237.799    34585.739
 66  TN+F          18381.270    198 37158.539    37575.492    37942.818
 67  TN+F+I        17444.925    199 35287.849    35711.253    36076.089
 68  TN+F+G4       16484.088    199 33366.175    33789.580    34154.416
 69  TN+F+I+G4     16434.375    200 33268.750    33698.697    34060.951
 70  TN+F+R2       16679.567    200 33759.134    34189.081    34551.335
 71  TN+F+R3       16494.503    202 33393.007    33836.315    34193.130
 72  TN+F+R4       16419.711    204 33247.421    33704.471    34055.467
 73  TN+F+R5       16402.297    206 33216.594    33687.776    34032.561
 74  TN+F+R6       16400.162    208 33216.324    33702.044    34040.213
 79  K3P           18432.266    195 37254.532    37652.657    38026.928
 80  K3P+I         17555.110    196 35502.219    35906.534    36278.577
 81  K3P+G4        16645.106    196 33682.212    34086.526    34458.569
 82  K3P+I+G4      16610.197    197 33614.393    34024.983    34394.711
 83  K3P+R2        16815.475    197 34024.949    34435.539    34805.267
 84  K3P+R3        16643.020    199 33684.040    34107.444    34472.280
 85  K3P+R4        16593.773    201 33589.547    34026.127    34385.709
 86  K3P+R5        16577.394    203 33560.787    34010.918    34364.871
 87  K3P+R6        16573.392    205 33556.783    34020.849    34368.789
 92  K3Pu+F        18330.508    198 37057.016    37473.968    37841.295
 93  K3Pu+F+I      17399.404    199 35196.808    35620.212    35985.048
 94  K3Pu+F+G4     16449.680    199 33297.359    33720.763    34085.599
 95  K3Pu+F+I+G4   16400.555    200 33201.110    33631.056    33993.311
 96  K3Pu+F+R2     16641.541    200 33683.082    34113.029    34475.283
 97  K3Pu+F+R3     16455.208    202 33314.416    33757.725    34114.540
 98  K3Pu+F+R4     16385.705    204 33179.410    33636.459    33987.455
 99  K3Pu+F+R5     16368.055    206 33148.109    33619.292    33964.076
100  K3Pu+F+R6     16366.249    208 33148.498    33634.219    33972.387
105  TPM2+F        18286.315    198 36968.631    37385.583    37752.910
106  TPM2+F+I      17358.379    199 35114.758    35538.162    35902.998
107  TPM2+F+G4     16422.806    199 33243.611    33667.015    34031.851
108  TPM2+F+I+G4   16376.561    200 33153.121    33583.068    33945.323
109  TPM2+F+R2     16616.254    200 33632.508    34062.454    34424.709
110  TPM2+F+R3     16427.765    202 33259.531    33702.839    34059.654
111  TPM2+F+R4     16362.701    204 33133.402    33590.451    33941.447
112  TPM2+F+R5     16346.215    206 33104.430    33575.612    33920.397
113  TPM2+F+R6     16344.329    208 33104.659    33590.379    33928.548
118  TPM2u+F       18286.315    198 36968.629    37385.582    37752.908
119  TPM2u+F+I     17358.378    199 35114.756    35538.160    35902.996
120  TPM2u+F+G4    16422.806    199 33243.612    33667.016    34031.852
121  TPM2u+F+I+G4  16376.574    200 33153.148    33583.094    33945.349
122  TPM2u+F+R2    16616.255    200 33632.509    34062.456    34424.710
123  TPM2u+F+R3    16427.584    202 33259.168    33702.476    34059.291
124  TPM2u+F+R4    16362.598    204 33133.197    33590.246    33941.242
125  TPM2u+F+R5    16346.214    206 33104.427    33575.610    33920.394
126  TPM2u+F+R6    16344.304    208 33104.607    33590.328    33928.497
131  TPM3+F        18328.432    198 37052.865    37469.817    37837.144
132  TPM3+F+I      17416.272    199 35230.543    35653.947    36018.783
133  TPM3+F+G4     16460.345    199 33318.689    33742.093    34106.929
134  TPM3+F+I+G4   16411.296    200 33222.592    33652.538    34014.793
135  TPM3+F+R2     16646.726    200 33693.452    34123.399    34485.653
136  TPM3+F+R3     16464.550    202 33333.101    33776.409    34133.224
137  TPM3+F+R4     16399.081    204 33206.161    33663.211    34014.207
138  TPM3+F+R5     16381.130    206 33174.259    33645.442    33990.226
139  TPM3+F+R6     16377.950    208 33171.899    33657.620    33995.788
144  TPM3u+F       18328.470    198 37052.941    37469.893    37837.220
145  TPM3u+F+I     17416.269    199 35230.538    35653.943    36018.778
146  TPM3u+F+G4    16460.344    199 33318.688    33742.092    34106.928
147  TPM3u+F+I+G4  16411.283    200 33222.566    33652.513    34014.767
148  TPM3u+F+R2    16646.695    200 33693.391    34123.337    34485.592
149  TPM3u+F+R3    16464.470    202 33332.941    33776.249    34133.064
150  TPM3u+F+R4    16399.020    204 33206.040    33663.090    34014.085
151  TPM3u+F+R5    16381.131    206 33174.261    33645.444    33990.229
152  TPM3u+F+R6    16377.905    208 33171.810    33657.530    33995.699
157  TIMe          18421.191    196 37234.382    37638.696    38010.739
158  TIMe+I        17540.899    197 35475.798    35886.388    36256.116
159  TIMe+G4       16635.698    197 33665.397    34075.986    34445.715
160  TIMe+I+G4     16601.851    198 33599.703    34016.655    34383.982
161  TIMe+R2       16806.316    198 34008.632    34425.584    34792.911
162  TIMe+R3       16634.097    200 33668.193    34098.140    34460.394
163  TIMe+R4       16584.031    202 33572.063    34015.371    34372.186
164  TIMe+R5       16566.625    204 33541.250    33998.299    34349.295
165  TIMe+R6       16564.195    206 33540.390    34011.572    34356.357
170  TIM+F         18327.164    199 37052.329    37475.733    37840.569
171  TIM+F+I       17394.532    200 35189.064    35619.010    35981.265
172  TIM+F+G4      16440.289    200 33280.578    33710.525    34072.779
173  TIM+F+I+G4    16393.607    201 33189.214    33625.794    33985.376
174  TIM+F+R2      16634.271    201 33670.542    34107.122    34466.704
175  TIM+F+R3      16446.826    203 33299.651    33749.782    34103.735
176  TIM+F+R4      16375.310    205 33160.619    33624.685    33972.625
177  TIM+F+R5      16357.617    207 33129.234    33607.634    33949.162
178  TIM+F+R6      16357.412    209 33132.825    33625.971    33960.675
183  TIM2e         18366.199    196 37124.399    37528.713    37900.756
184  TIM2e+I       17502.561    197 35399.123    35809.712    36179.441
185  TIM2e+G4      16619.998    197 33633.997    34044.586    34414.315
186  TIM2e+I+G4    16592.536    198 33581.071    33998.024    34365.350
187  TIM2e+R2      16789.695    198 33975.389    34392.342    34759.668
188  TIM2e+R3      16618.860    200 33637.720    34067.667    34429.921
189  TIM2e+R4      16576.250    202 33556.500    33999.808    34356.623
190  TIM2e+R5      16554.976    204 33517.952    33975.001    34325.997
191  TIM2e+R6      16554.092    206 33520.185    33991.367    34336.152
196  TIM2+F        18283.988    199 36965.975    37389.379    37754.215
197  TIM2+F+I      17354.905    200 35109.809    35539.756    35902.010
198  TIM2+F+G4     16412.903    200 33225.805    33655.752    34018.006
199  TIM2+F+I+G4   16372.675    201 33147.350    33583.931    33943.512
200  TIM2+F+R2     16609.552    201 33621.104    34057.685    34417.266
201  TIM2+F+R3     16420.020    203 33246.039    33696.170    34050.123
202  TIM2+F+R4     16354.365    205 33118.730    33582.796    33930.736
203  TIM2+F+R5     16332.887    207 33079.774    33558.174    33899.702
204  TIM2+F+R6     16332.783    209 33083.566    33576.712    33911.416
209  TIM3e         18306.803    196 37005.605    37409.919    37781.962
210  TIM3e+I       17474.723    197 35343.445    35754.035    36123.763
211  TIM3e+G4      16555.480    197 33504.960    33915.550    34285.278
212  TIM3e+I+G4    16513.812    198 33423.624    33840.577    34207.903
213  TIM3e+R2      16718.707    198 33833.415    34250.367    34617.694
214  TIM3e+R3      16555.100    200 33510.200    33940.147    34302.401
215  TIM3e+R4      16500.831    202 33405.662    33848.970    34205.785
216  TIM3e+R5      16482.693    204 33373.386    33830.435    34181.431
217  TIM3e+R6      16482.650    206 33377.300    33848.482    34193.267
222  TIM3+F        18324.419    199 37046.838    37470.242    37835.078
223  TIM3+F+I      17411.087    200 35222.173    35652.120    36014.374
224  TIM3+F+G4     16454.013    200 33308.026    33737.972    34100.227
225  TIM3+F+I+G4   16404.895    201 33211.789    33648.370    34007.951
226  TIM3+F+R2     16641.067    201 33684.134    34120.715    34480.296
227  TIM3+F+R3     16457.422    203 33320.845    33770.975    34124.929
228  TIM3+F+R4     16390.467    205 33190.935    33655.001    34002.941
229  TIM3+F+R5     16368.323    207 33150.645    33629.045    33970.573
230  TIM3+F+R6     16368.193    209 33154.386    33647.533    33982.237
235  TVMe          18153.293    197 36700.585    37111.175    37480.903
236  TVMe+I        17330.193    198 35056.387    35473.339    35840.666
237  TVMe+G4       16402.072    198 33200.143    33617.095    33984.422
238  TVMe+I+G4     16370.952    199 33139.904    33563.309    33928.144
239  TVMe+R2       16570.496    199 33538.992    33962.396    34327.232
240  TVMe+R3       16403.304    201 33208.607    33645.188    34004.769
241  TVMe+R4       16356.376    203 33118.753    33568.883    33922.837
242  TVMe+R5       16340.383    205 33090.767    33554.833    33902.773
243  TVMe+R6       16340.304    207 33094.608    33573.008    33914.536
248  TVM+F         18230.263    200 36860.525    37290.472    37652.726
249  TVM+F+I       17324.262    201 35050.524    35487.105    35846.686
250  TVM+F+G4      16386.932    201 33175.865    33612.446    33972.027
251  TVM+F+I+G4    16344.859    202 33093.718    33537.026    33893.841
252  TVM+F+R2      16574.379    202 33552.759    33996.067    34352.882
253  TVM+F+R3      16388.958    204 33185.916    33642.965    33993.961
254  TVM+F+R4      16330.994    206 33073.989    33545.171    33889.956
255  TVM+F+R5      16310.964    208 33037.929    33523.649    33861.818
256  TVM+F+R6      16310.941    210 33041.883    33542.561    33873.694
261  SYM           18142.532    198 36681.064    37098.016    37465.343
262  SYM+I         17315.569    199 35029.138    35452.543    35817.378
263  SYM+G4        16398.344    199 33194.688    33618.092    33982.928
264  SYM+I+G4      16366.677    200 33133.355    33563.302    33925.556
265  SYM+R2        16564.826    200 33529.652    33959.598    34321.853
266  SYM+R3        16398.573    202 33201.146    33644.454    34001.269
267  SYM+R4        16353.302    204 33114.603    33571.652    33922.648
268  SYM+R5        16336.728    206 33085.456    33556.638    33901.423
269  SYM+R6        16336.520    208 33089.040    33574.760    33912.929
274  GTR+F         18227.300    201 36856.600    37293.181    37652.762
275  GTR+F+I       17320.281    202 35044.562    35487.870    35844.685
276  GTR+F+G4      16379.610    202 33163.220    33606.528    33963.343
277  GTR+F+I+G4    16339.379    203 33084.758    33534.888    33888.842
278  GTR+F+R2      16568.444    203 33542.888    33993.018    34346.972
279  GTR+F+R3      16382.064    205 33174.128    33638.194    33986.134
280  GTR+F+R4      16322.959    207 33059.919    33538.319    33879.847
281  GTR+F+R5      16302.328    209 33022.655    33515.801    33850.505
282  GTR+F+R6      16302.037    211 33026.074    33534.392    33861.846
Akaike Information Criterion:           GTR+F+R5
Corrected Akaike Information Criterion: GTR+F+R5
Bayesian Information Criterion:         GTR+F+R5
Best-fit model: GTR+F+R5 chosen according to BIC

All model information printed to 16S-aligned-mafft.fasta.model.gz
CPU time for ModelFinder: 60.425 seconds (0h:1m:0s)
Wall-clock time for ModelFinder: 60.546 seconds (0h:1m:0s)

NOTE: 5 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -17058.435
2. Current log-likelihood: -16372.749
3. Current log-likelihood: -16356.384
4. Current log-likelihood: -16350.509
5. Current log-likelihood: -16345.599
6. Current log-likelihood: -16341.342
7. Current log-likelihood: -16337.639
8. Current log-likelihood: -16334.167
9. Current log-likelihood: -16331.978
10. Current log-likelihood: -16330.060
11. Current log-likelihood: -16328.922
12. Current log-likelihood: -16327.643
13. Current log-likelihood: -16326.960
14. Current log-likelihood: -16326.069
15. Current log-likelihood: -16325.590
16. Current log-likelihood: -16324.867
17. Current log-likelihood: -16324.647
18. Current log-likelihood: -16324.388
19. Current log-likelihood: -16324.080
20. Current log-likelihood: -16323.886
21. Current log-likelihood: -16323.649
22. Current log-likelihood: -16323.529
23. Current log-likelihood: -16323.303
Optimal log-likelihood: -16323.195
Rate parameters:  A-C: 1.23115  A-G: 5.10286  A-T: 2.66445  C-G: 0.59213  C-T: 6.39806  G-T: 1.00000
Base frequencies:  A: 0.313  C: 0.142  G: 0.228  T: 0.317
Site proportion and rates:  (0.338,0.064) (0.154,0.372) (0.130,0.884) (0.298,1.780) (0.080,3.447)
Parameters optimization took 23 rounds (2.787 sec)
Computing ML distances based on estimated model parameters... 0.066 sec
Computing BIONJ tree...
0.006 seconds
Log-likelihood of BIONJ tree: -16271.210
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 0.584 second
Computing log-likelihood of 98 initial trees ... 1.532 seconds
Current best score: -16247.193

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -16195.413
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 3: -16189.130
Iteration 10 / LogL: -16194.866 / Time: 0h:0m:9s
Iteration 20 / LogL: -16209.725 / Time: 0h:0m:11s
Finish initializing candidate tree set (20)
Current best tree score: -16189.130 / CPU time: 8.241
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 22: -16185.393
Estimate model parameters (epsilon = 0.100)
UPDATE BEST LOG-LIKELIHOOD: -16184.610
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 29: -16183.401
Iteration 30 / LogL: -16184.099 / Time: 0h:0m:14s (0h:0m:51s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 40: -16182.903
Iteration 40 / LogL: -16182.903 / Time: 0h:0m:17s (0h:0m:43s left)
Iteration 50 / LogL: -16199.570 / Time: 0h:0m:18s (0h:0m:34s left)
Iteration 60 / LogL: -16190.700 / Time: 0h:0m:21s (0h:0m:28s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 63: -16182.541
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 69: -16182.211
Iteration 70 / LogL: -16182.682 / Time: 0h:0m:23s (0h:0m:33s left)
UPDATE BEST LOG-LIKELIHOOD: -16182.203
Iteration 80 / LogL: -16185.493 / Time: 0h:0m:25s (0h:0m:28s left)
Iteration 90 / LogL: -16214.157 / Time: 0h:0m:27s (0h:0m:24s left)
Iteration 100 / LogL: -16182.216 / Time: 0h:0m:29s (0h:0m:20s left)
Iteration 110 / LogL: -16191.187 / Time: 0h:0m:31s (0h:0m:16s left)
UPDATE BEST LOG-LIKELIHOOD: -16182.201
UPDATE BEST LOG-LIKELIHOOD: -16182.201
Iteration 120 / LogL: -16191.223 / Time: 0h:0m:33s (0h:0m:13s left)
Iteration 130 / LogL: -16182.330 / Time: 0h:0m:34s (0h:0m:10s left)
UPDATE BEST LOG-LIKELIHOOD: -16182.199
Iteration 140 / LogL: -16194.207 / Time: 0h:0m:36s (0h:0m:7s left)
Iteration 150 / LogL: -16184.892 / Time: 0h:0m:38s (0h:0m:4s left)
Iteration 160 / LogL: -16182.356 / Time: 0h:0m:40s (0h:0m:2s left)
Iteration 170 / LogL: -16184.650 / Time: 0h:0m:42s (0h:0m:-1s left)
TREE SEARCH COMPLETED AFTER 170 ITERATIONS / Time: 0h:0m:42s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -16182.199
2. Current log-likelihood: -16182.113
3. Current log-likelihood: -16182.052
4. Current log-likelihood: -16181.967
5. Current log-likelihood: -16181.916
6. Current log-likelihood: -16181.838
7. Current log-likelihood: -16181.794
8. Current log-likelihood: -16181.726
9. Current log-likelihood: -16181.696
10. Current log-likelihood: -16181.625
11. Current log-likelihood: -16181.609
12. Current log-likelihood: -16181.560
13. Current log-likelihood: -16181.548
14. Current log-likelihood: -16181.511
15. Current log-likelihood: -16181.499
16. Current log-likelihood: -16181.464
17. Current log-likelihood: -16181.452
18. Current log-likelihood: -16181.420
Optimal log-likelihood: -16181.404
Rate parameters:  A-C: 1.40801  A-G: 6.11492  A-T: 3.11074  C-G: 0.64074  C-T: 7.59328  G-T: 1.00000
Base frequencies:  A: 0.313  C: 0.142  G: 0.228  T: 0.317
Site proportion and rates:  (0.287,0.034) (0.199,0.272) (0.165,0.892) (0.287,1.845) (0.062,4.164)
Parameters optimization took 18 rounds (1.878 sec)
BEST SCORE FOUND : -16181.404
Total tree length: 17.802

Total number of iterations: 170
CPU time used for tree search: 39.824 sec (0h:0m:39s)
Wall-clock time used for tree search: 39.896 sec (0h:0m:39s)
Total CPU time used: 44.610 sec (0h:0m:44s)
Total wall-clock time used: 44.695 sec (0h:0m:44s)

Analysis results written to: 
  IQ-TREE report:                16S-aligned-mafft.fasta.iqtree
  Maximum-likelihood tree:       16S-aligned-mafft.fasta.treefile
  Likelihood distances:          16S-aligned-mafft.fasta.mldist
  Screen log file:               16S-aligned-mafft.fasta.log
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/16S-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies/Gene_Trees

#### 18S Gene Tree
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s 18S-aligned-mafft.fasta>
-running 18S Mafft alignment
Output: 
Reading alignment file 18S-aligned-mafft.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
WARNING: 1 sites contain only gaps or ambiguous characters.
Alignment has 191 sequences with 1762 columns, 776 distinct patterns
269 parsimony-informative, 159 singleton sites, 1334 constant sites
                         Gap/Ambiguity  Composition  p-value
   1  Pseudoepedanus            59.65%    passed     46.09%
   2  Agoristenidae_DNA105839    2.50%    passed     96.46%
   3  Trinella                   0.62%    passed     94.94%
   4  Arulla_DNA102666           0.74%    passed     93.88%
   5  Assamiidae_DNA104069       0.79%    passed     97.01%
   6  Assamiidae_DNA104857       0.62%    passed     93.93%
   7  Paktongius                 0.57%    passed     94.62%
   8  Chilon                    58.00%    passed     80.62%
   9  Paraselenca                0.57%    passed     98.13%
  10  Assamiidae_DNA104859       1.08%    passed     93.33%
  11  Assamiidae_DNA104858       0.68%    passed     98.22%
  12  Baculigerus                0.57%    passed     98.03%
  13  Fijicolana                 0.57%    passed     97.50%
  14  Stygnomma                 62.66%    passed     64.71%
  15  Pellobunus                 0.57%    passed     97.52%
  16  Icaleptes_DNA104056        0.57%    passed     97.50%
  17  Icaleptes_DNA104053        0.57%    passed     96.71%
  18  Guasinia                   0.57%    passed     97.84%
  19  Metabiantes_DNA100703      0.57%    passed     95.43%
  20  Metabiantes_DNA100704      0.57%    passed     95.82%
  21  Stenostygnus_DNA104847     0.68%    passed     97.77%
  22  Stenostygnus_DNA104849     0.57%    passed     95.51%
  23  Stenostygnus_DNA104850     0.85%    passed     96.94%
  24  Stenostygnus_DNA104848     0.57%    passed     97.24%
  25  Ethobunus                  0.57%    passed     95.73%
  26  Fissiphallius_DNA104055    0.57%    passed     96.26%
  27  Fissiphallius_DNA104057    5.11%    passed     88.48%
  28  Zalmoxis                   0.68%    passed     95.49%
  29  Escadabius                 0.57%    passed     98.15%
  30  Metabiantes_DNA100335      0.57%    passed     93.98%
  31  Biantidae_DNA105668        0.57%    passed     94.78%
  32  Samoidae                   0.57%    passed     97.32%
  33  Caenoncopus                0.57%    passed     96.15%
  34  Gnomulus                   1.19%    passed     94.20%
  35  Palaeoncopus               0.57%    passed     95.01%
  36  Sandokan                   0.57%    passed     93.19%
  37  Conomma                    0.57%    passed     97.80%
  38  Mitraceras                 0.57%    passed     97.34%
  39  Haasus_sp                  0.57%    passed     97.80%
  40  Haasus_judaeus             4.26%    passed     92.48%
  41  Jarmilana                  2.27%    passed     93.64%
  42  Maiorerus                  0.62%    passed     94.55%
  43  Bishopella                 0.57%    passed     97.62%
  44  Epedanidae_DNA104861       0.57%    passed     96.78%
  45  Epedanidae_DNA104862       0.79%    passed     94.36%
  46  Cynortula                  2.04%    passed     90.20%
  47  Glysterus                  0.57%    passed     95.77%
  48  Goniosoma                  0.62%    passed     94.67%
  49  Megapachylus               8.34%    passed     66.26%
  50  Heterocranaus              0.57%    passed     94.42%
  51  Santinezia                 0.57%    passed     93.04%
  52  Zygopachylus               0.57%    passed     94.48%
  53  Metalibitia                0.57%    passed     86.87%
  54  Martensiellus              0.62%    passed     97.17%
  55  Santobius_DNA104930        0.91%    passed     97.09%
  56  Santobius_DNA104931        0.57%    passed     98.12%
  57  Dongmoa                    0.62%    passed     91.32%
  58  Lomanius_DNA104935         0.62%    passed     94.09%
  59  Hoplobunus                 0.57%    passed     95.77%
  60  Stygnopsidae_DNA103882     0.57%    passed     97.42%
  61  Stygnopsidae_DNA104855     2.44%    passed     97.07%
  62  Stygnopsidae_DNA104856     0.57%    passed     96.15%
  63  Minuella                   5.11%    passed     96.79%
  64  Trionyxella                5.11%    passed     92.30%
  65  Karos                      0.57%    passed     95.81%
  66  Equitius                   0.74%    passed     97.02%
  67  Larifuga                  26.73%    passed     92.34%
  68  Triaenonychidae            0.57%    passed     97.92%
  69  Triaenobunus               0.57%    passed     97.34%
  70  Rostromontia               0.57%    passed     99.02%
  71  IC_DNA102668               0.79%    passed     96.49%
  72  IC_DNA102669               0.68%    passed     96.92%
  73  Zalmoxida                  0.74%    passed     98.72%
  74  IC_DNA103729               1.14%    passed     97.59%
  75  IC_DNA104070               0.79%    passed     98.80%
  76  IC_DNA104071               0.74%    passed     98.45%
  77  Synthetonychia            24.52%    passed     82.53%
  78  IC_DNA103572              21.85%    passed     97.43%
  79  Bunofagea                  4.54%    passed     75.32%
  80  Scotolemon                 1.31%    passed     90.79%
  81  Scotolemon_lespesi         0.74%    passed     97.32%
  82  Vonones                   46.88%    passed     80.41%
  83  Erebomaster                0.62%    passed     91.91%
  84  Holoscotolemon             0.57%    passed     96.07%
  85  Peltonychia                0.57%    passed     96.79%
  86  Trojanella                 0.57%    passed     95.82%
  87  Epedanidae_DNA104062       0.62%    passed     93.58%
  88  Epedanidae_DNA104068       0.85%    passed     92.58%
  89  Epedanidae_DNA104066       0.74%    passed     92.52%
  90  Tithaeus                   0.62%    passed     94.03%
  91  Op106_Toccolus             0.17%    passed     94.08%
  92  Op107_Nanepedanus_rufus    0.17%    passed     95.97%
  93  Op104_Dibuninae            0.17%    passed     96.28%
  94  Op105_Dibunus              0.17%    passed     96.60%
  95  As020_Gab                 51.19%    passed     69.04%
  96  As026_Gab                 56.41%    passed     83.20%
  97  As041_Gab                 49.04%    passed     83.53%
  98  As030_Gab                 55.96%    passed     89.06%
  99  As038_Gab                 49.32%    passed     89.61%
 100  As022_Gab                 54.37%    passed     62.63%
 101  As098_Phil                50.68%    passed     88.23%
 102  As105_Phil                49.77%    passed     88.82%
 103  As011_Gab                 59.42%    passed     93.71%
 104  As025_Gab                 60.39%    passed     93.13%
 105  As070_Lib                 49.04%    passed     87.14%
 106  As108_Laos                49.04%    passed     85.77%
 107  As028_Gab                 55.39%    passed     72.51%
 108  As114_Thai                49.15%    passed     82.88%
 109  As111_Phil                49.66%    passed     84.60%
 110  As083_Cam                 49.21%    passed     86.51%
 111  As084_Cam                 50.00%    passed     88.09%
 112  As085_Cam                 49.26%    passed     86.71%
 113  As109_Thai                49.49%    passed     86.78%
 114  As121_Laos                49.09%    passed     85.86%
 115  As125_Thai                49.32%    passed     86.42%
 116  As132_Mala                49.38%    passed     81.27%
 117  As039_Gab                 49.32%    passed     82.84%
 118  As027_Gab                 58.85%    passed     78.65%
 119  As058_Gab                 49.21%    passed     89.37%
 120  As085_WAus                49.15%    passed     82.07%
 121  As092_EAus                49.32%    passed     82.84%
 122  Op103_Bupares              0.17%    passed     97.66%
 123  As120_Thai                49.32%    passed     89.49%
 124  As012_Gab                 62.03%    passed     59.96%
 125  As097_PNG                 49.09%    passed     91.34%
 126  As101_NAus                50.06%    passed     97.41%
 127  As133_Indo                49.15%    passed     87.79%
 128  As034_Cam                 50.74%    passed     94.68%
 129  As040_Gab                 49.43%    passed     83.07%
 130  As059_Gab                 49.09%    passed     88.17%
 131  As029_Gab                 51.36%    passed     86.28%
 132  As057_Gab                 49.32%    passed     85.89%
 133  As061_Lib                 49.38%    passed     85.42%
 134  As094_PNG                 49.60%    passed     82.96%
 135  As096_PNG                 49.15%    passed     89.04%
 136  As130_Viet                49.09%    passed     86.86%
 137  As117_Thai                49.38%    passed     85.96%
 138  As116_Thai                49.21%    passed     87.55%
 139  As124_Indo                49.32%    passed     84.81%
 140  As134_Viet                49.43%    passed     86.26%
 141  As089_EAus                49.38%    passed     86.89%
 142  As090_EAus                49.32%    passed     84.86%
 143  As103_NAus                51.08%    passed     89.38%
 144  As045_Cam                 49.43%    passed     84.80%
 145  As069_Lib                 49.72%    passed     84.08%
 146  As021_Gab                 60.61%    passed     90.03%
 147  As127_Laos                49.43%    passed     93.57%
 148  As087_WAus                49.38%    passed     84.02%
 149  As032_Gab                 52.16%    passed     87.39%
 150  As056_Gab                 49.49%    passed     82.96%
 151  As081_Indo                49.32%    passed     89.49%
 152  As095_PNG                 49.77%    passed     92.22%
 153  As082_Cam                 49.49%    passed     84.48%
 154  As080_Indo                49.72%    passed     87.38%
 155  As072_Lib                 49.55%    passed     89.63%
 156  As110_Laos                49.66%    passed     78.52%
 157  As009_Cam                 54.09%    passed     82.08%
 158  As010_Gab                 57.78%    passed     88.03%
 159  As016_Gab                 59.70%    passed     92.23%
 160  As018_Gab                 56.07%    passed     80.42%
 161  As050_Cam                 49.15%    passed     86.81%
 162  As031_Gab                 52.95%    passed     59.36%
 163  As102_NAus                50.17%    passed     87.30%
 164  As017_Gab                 52.27%    passed     79.72%
 165  Op049_Beloniscus           1.70%    passed     97.70%
 166  As135_Viet                49.15%    passed     84.50%
 167  As106_Thai                49.21%    passed     82.13%
 168  As115_Thai                49.55%    passed     92.84%
 169  As118_Thai                49.21%    passed     83.14%
 170  As122_Laos                61.29%    passed     82.46%
 171  As129_Laos                49.60%    passed     88.43%
 172  As136_Viet                49.15%    passed     91.24%
 173  As140_Viet                49.09%    passed     90.95%
 174  As123_Viet                49.38%    passed     88.50%
 175  As126_Laos                49.49%    passed     91.41%
 176  As119_Indo                49.09%    passed     84.73%
 177  As131_Mala                50.85%    passed     79.18%
 178  Kimula                    28.15%    passed     79.07%
 179  Fumontana                  0.57%    passed     65.37%
 180  As104_EAus                51.93%    passed     87.17%
 181  Zuma                       0.62%    passed     95.00%
 182  Stygnoplus                 0.68%    passed     95.73%
 183  Caddo                      0.79%    passed     99.61%
 184  Pantopsalis                0.74%    passed     96.88%
 185  Protolophus                5.96%    passed     96.72%
 186  Theromaster                0.62%    passed     82.41%
 187  As099_Dibunus             48.52%    passed     88.03%
 188  Troglosiro                 0.79%    passed     98.29%
 189  Dendrolasma                0.79%    passed     88.40%
 190  Trogulus                   0.74%    passed     80.55%
 191  Hesperonemastoma           0.62%    passed     95.76%
WARNING: 29 sequences contain more than 50% gaps/ambiguity
****  TOTAL                     24.47%  0 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.014 seconds
NOTE: ModelFinder requires 48 MB RAM!
ModelFinder will test 286 DNA models (sample size: 1762) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  JC            11405.676    379 23569.352    23777.775    25644.076
  2  JC+I          10527.104    380 21814.209    22023.883    23894.407
  3  JC+G4         10390.169    380 21540.337    21750.012    23620.535
  4  JC+I+G4       10258.178    381 21278.357    21489.287    23364.029
  5  JC+R2         10273.165    381 21308.330    21519.261    23394.002
  6  JC+R3         10182.056    383 21130.112    21343.569    23226.733
  7  JC+R4         10181.603    385 21133.207    21349.210    23240.776
 14  F81+F         11417.280    382 23598.560    23810.751    25689.706
 15  F81+F+I       10541.354    383 21848.708    22062.165    23945.329
 16  F81+F+G4      10404.408    383 21574.815    21788.273    23671.436
 17  F81+F+I+G4    10273.849    384 21315.699    21530.427    23417.794
 18  F81+F+R2      10288.128    384 21344.255    21558.983    23446.350
 19  F81+F+R3      10197.766    386 21167.532    21384.815    23280.575
 20  F81+F+R4      10197.454    388 21170.908    21390.765    23294.900
 27  K2P           11160.727    380 23081.455    23291.129    25161.652
 28  K2P+I         10277.324    381 21316.649    21527.579    23402.321
 29  K2P+G4        10135.944    381 21033.888    21244.819    23119.560
 30  K2P+I+G4      10001.658    382 20767.315    20979.507    22858.461
 31  K2P+R2        10017.514    382 20799.029    21011.220    22890.175
 32  K2P+R3        9922.906     384 20613.812    20828.539    22715.906
 33  K2P+R4        9922.520     386 20617.041    20834.324    22730.084
 40  HKY+F         11173.975    383 23113.950    23327.407    25210.571
 41  HKY+F+I       10288.899    384 21345.798    21560.526    23447.893
 42  HKY+F+G4      10145.691    384 21059.382    21274.110    23161.477
 43  HKY+F+I+G4    10011.767    385 20793.535    21009.538    22901.103
 44  HKY+F+R2      10027.638    385 20825.275    21041.278    22932.844
 45  HKY+F+R3      9932.052     387 20638.105    20856.673    22756.622
 46  HKY+F+R4      9931.720     389 20641.439    20862.591    22770.905
 53  TNe           11112.502    381 22987.004    23197.935    25072.676
 54  TNe+I         10248.450    382 21260.899    21473.091    23352.046
 55  TNe+G4        10095.917    382 20955.834    21168.026    23046.980
 56  TNe+I+G4      9973.158     383 20712.315    20925.772    22808.936
 57  TNe+R2        9989.822     383 20745.644    20959.101    22842.264
 58  TNe+R3        9900.189     385 20570.378    20786.381    22677.947
 59  TNe+R4        9899.996     387 20573.992    20792.560    22692.509
 66  TN+F          11108.037    384 22984.074    23198.802    25086.169
 67  TN+F+I        10246.868    385 21263.736    21479.739    23371.305
 68  TN+F+G4       10088.433    385 20946.866    21162.869    23054.435
 69  TN+F+I+G4     9968.414     386 20708.829    20926.111    22821.872
 70  TN+F+R2       9985.938     386 20743.876    20961.159    22856.919
 71  TN+F+R3       9896.485     388 20568.970    20788.827    22692.962
 72  TN+F+R4       9896.330     390 20572.660    20795.111    22707.600
 79  K3P           11158.375    381 23078.750    23289.681    25164.423
 80  K3P+I         10274.764    382 21313.527    21525.719    23404.673
 81  K3P+G4        10133.103    382 21030.205    21242.397    23121.351
 82  K3P+I+G4      9998.972     383 20763.944    20977.401    22860.564
 83  K3P+R2        10014.619    383 20795.238    21008.695    22891.858
 84  K3P+R3        9920.134     385 20610.268    20826.271    22717.837
 85  K3P+R4        9919.779     387 20613.557    20832.125    22732.075
 92  K3Pu+F        11171.695    384 23111.391    23326.119    25213.486
 93  K3Pu+F+I      10286.395    385 21342.791    21558.793    23450.359
 94  K3Pu+F+G4     10143.006    385 21056.012    21272.015    23163.581
 95  K3Pu+F+I+G4   10009.176    386 20790.353    21007.636    22903.396
 96  K3Pu+F+R2     10024.866    386 20821.732    21039.014    22934.775
 97  K3Pu+F+R3     9929.473     388 20634.946    20854.803    22758.938
 98  K3Pu+F+R4     9929.112     390 20638.224    20860.675    22773.164
105  TPM2+F        11172.830    384 23113.659    23328.387    25215.754
106  TPM2+F+I      10286.756    385 21343.512    21559.514    23451.080
107  TPM2+F+G4     10143.265    385 21056.529    21272.532    23164.098
108  TPM2+F+I+G4   10009.384    386 20790.769    21008.051    22903.812
109  TPM2+F+R2     10026.080    386 20824.159    21041.442    22937.202
110  TPM2+F+R3     9929.486     388 20634.972    20854.829    22758.963
111  TPM2+F+R4     9929.140     390 20638.280    20860.731    22773.220
118  TPM2u+F       11172.828    384 23113.655    23328.383    25215.750
119  TPM2u+F+I     10286.754    385 21343.508    21559.511    23451.077
120  TPM2u+F+G4    10143.264    385 21056.527    21272.530    23164.096
121  TPM2u+F+I+G4  10009.389    386 20790.778    21008.061    22903.821
122  TPM2u+F+R2    10026.075    386 20824.150    21041.433    22937.193
123  TPM2u+F+R3    9929.484     388 20634.968    20854.825    22758.959
124  TPM2u+F+R4    9929.134     390 20638.268    20860.719    22773.208
131  TPM3+F        11171.948    384 23111.896    23326.624    25213.991
132  TPM3+F+I      10285.041    385 21340.082    21556.085    23447.651
133  TPM3+F+G4     10144.038    385 21058.075    21274.078    23165.644
134  TPM3+F+I+G4   10008.762    386 20789.523    21006.806    22902.566
135  TPM3+F+R2     10026.023    386 20824.045    21041.328    22937.088
136  TPM3+F+R3     9930.179     388 20636.358    20856.215    22760.350
137  TPM3+F+R4     9929.804     390 20639.607    20862.058    22774.547
144  TPM3u+F       11171.950    384 23111.900    23326.628    25213.995
145  TPM3u+F+I     10285.043    385 21340.086    21556.089    23447.655
146  TPM3u+F+G4    10144.038    385 21058.077    21274.080    23165.646
147  TPM3u+F+I+G4  10008.761    386 20789.522    21006.805    22902.565
148  TPM3u+F+R2    10026.028    386 20824.057    21041.340    22937.100
149  TPM3u+F+R3    9930.168     388 20636.337    20856.194    22760.328
150  TPM3u+F+R4    9929.798     390 20639.595    20862.046    22774.535
157  TIMe          11110.127    382 22984.253    23196.445    25075.399
158  TIMe+I        10245.884    383 21257.768    21471.225    23354.389
159  TIMe+G4       10093.107    383 20952.213    21165.670    23048.834
160  TIMe+I+G4     9970.442     384 20708.883    20923.611    22810.978
161  TIMe+R2       9986.958     384 20741.915    20956.643    22844.010
162  TIMe+R3       9897.484     386 20566.968    20784.250    22680.011
163  TIMe+R4       9897.244     388 20570.488    20790.346    22694.480
170  TIM+F         11105.756    385 22981.511    23197.514    25089.080
171  TIM+F+I       10244.380    386 21260.761    21478.044    23373.804
172  TIM+F+G4      10085.819    386 20943.639    21160.922    23056.682
173  TIM+F+I+G4    9965.839     387 20705.678    20924.246    22824.196
174  TIM+F+R2      9983.227     387 20740.453    20959.021    22858.970
175  TIM+F+R3      9893.942     389 20565.884    20787.036    22695.350
176  TIM+F+R4      9893.736     391 20569.471    20793.226    22709.885
183  TIM2e         11111.912    382 22987.823    23200.015    25078.969
184  TIM2e+I       10247.118    383 21260.236    21473.694    23356.857
185  TIM2e+G4      10093.897    383 20953.794    21167.251    23050.414
186  TIM2e+I+G4    9971.507     384 20711.015    20925.742    22813.109
187  TIM2e+R2      9988.703     384 20745.405    20960.133    22847.500
188  TIM2e+R3      9898.852     386 20569.704    20786.987    22682.747
189  TIM2e+R4      9898.598     388 20573.196    20793.054    22697.188
196  TIM2+F        11106.895    385 22983.790    23199.792    25091.358
197  TIM2+F+I      10244.677    386 21261.353    21478.636    23374.396
198  TIM2+F+G4     10085.146    386 20942.292    21159.575    23055.335
199  TIM2+F+I+G4   9965.575     387 20705.150    20923.718    22823.667
200  TIM2+F+R2     9983.733     387 20741.466    20960.034    22859.983
201  TIM2+F+R3     9894.101     389 20566.202    20787.353    22695.667
202  TIM2+F+R4     9893.893     391 20569.786    20793.541    22710.200
209  TIM3e         11109.354    382 22982.709    23194.900    25073.855
210  TIM3e+I       10243.304    383 21252.608    21466.065    23349.228
211  TIM3e+G4      10092.769    383 20951.538    21164.995    23048.158
212  TIM3e+I+G4    9968.858     384 20705.715    20920.443    22807.810
213  TIM3e+R2      9986.650     384 20741.300    20956.028    22843.395
214  TIM3e+R3      9897.455     386 20566.911    20784.194    22679.954
215  TIM3e+R4      9897.136     388 20570.272    20790.129    22694.263
222  TIM3+F        11106.035    385 22982.069    23198.072    25089.638
223  TIM3+F+I      10243.088    386 21258.177    21475.459    23371.220
224  TIM3+F+G4     10086.904    386 20945.808    21163.091    23058.851
225  TIM3+F+I+G4   9966.039     387 20706.079    20924.647    22824.596
226  TIM3+F+R2     9984.197     387 20742.395    20960.962    22860.912
227  TIM3+F+R3     9895.087     389 20568.174    20789.325    22697.639
228  TIM3+F+R4     9894.826     391 20571.652    20795.406    22712.066
235  TVMe          11154.207    383 23074.413    23287.870    25171.034
236  TVMe+I        10267.932    384 21303.864    21518.592    23405.959
237  TVMe+G4       10128.489    384 21024.978    21239.706    23127.073
238  TVMe+I+G4     9992.728     385 20755.457    20971.460    22863.026
239  TVMe+R2       10010.772    385 20791.544    21007.547    22899.113
240  TVMe+R3       9915.464     387 20604.929    20823.496    22723.446
241  TVMe+R4       9915.002     389 20608.004    20829.156    22737.470
248  TVM+F         11168.155    386 23108.310    23325.593    25221.353
249  TVM+F+I       10280.127    387 21334.254    21552.822    23452.771
250  TVM+F+G4      10138.600    387 21051.201    21269.768    23169.718
251  TVM+F+I+G4    10003.415    388 20782.830    21002.688    22906.822
252  TVM+F+R2      10021.394    388 20818.787    21038.645    22942.779
253  TVM+F+R3      9924.875     390 20629.749    20852.200    22764.689
254  TVM+F+R4      9924.414     392 20632.827    20857.891    22778.716
261  SYM           11105.967    384 22979.933    23194.661    25082.028
262  SYM+I         10239.098    385 21248.195    21464.198    23355.764
263  SYM+G4        10087.620    385 20945.240    21161.243    23052.809
264  SYM+I+G4      9964.218     386 20700.436    20917.719    22813.479
265  SYM+R2        9982.301     386 20736.602    20953.885    22849.645
266  SYM+R3        9893.046     388 20562.092    20781.949    22686.084
267  SYM+R4        9892.687     390 20565.374    20787.825    22700.314
274  GTR+F         11102.224    387 22978.448    23197.015    25096.965
275  GTR+F+I       10238.134    388 21252.268    21472.125    23376.259
276  GTR+F+G4      10080.595    388 20937.189    21157.047    23061.181
277  GTR+F+I+G4    9960.267     389 20698.534    20919.685    22828.000
278  GTR+F+R2      9978.885     389 20735.771    20956.923    22865.237
279  GTR+F+R3      9889.696     391 20561.392    20785.146    22701.806
280  GTR+F+R4      9889.410     393 20564.820    20791.197    22716.183
Akaike Information Criterion:           GTR+F+R3
Corrected Akaike Information Criterion: SYM+R3
Bayesian Information Criterion:         TNe+R3
Best-fit model: TNe+R3 chosen according to BIC

All model information printed to 18S-aligned-mafft.fasta.model.gz
CPU time for ModelFinder: 57.277 seconds (0h:0m:57s)
Wall-clock time for ModelFinder: 57.437 seconds (0h:0m:57s)

NOTE: 14 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -10861.021
2. Current log-likelihood: -10013.523
3. Current log-likelihood: -9974.484
4. Current log-likelihood: -9965.499
5. Current log-likelihood: -9954.484
6. Current log-likelihood: -9939.786
7. Current log-likelihood: -9922.545
8. Current log-likelihood: -9908.998
9. Current log-likelihood: -9903.751
10. Current log-likelihood: -9902.331
11. Current log-likelihood: -9901.637
12. Current log-likelihood: -9901.215
13. Current log-likelihood: -9900.918
14. Current log-likelihood: -9900.769
15. Current log-likelihood: -9900.647
Optimal log-likelihood: -9900.500
Rate parameters:  A-C: 1.00000  A-G: 2.71558  A-T: 1.00000  C-G: 1.00000  C-T: 5.15208  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Site proportion and rates:  (0.854,0.172) (0.114,3.344) (0.032,14.761)
Parameters optimization took 15 rounds (3.328 sec)
Computing ML distances based on estimated model parameters... 0.305 sec
WARNING: Some pairwise ML distances are too long (saturated)
Computing BIONJ tree...
0.016 seconds
Log-likelihood of BIONJ tree: -10182.875
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 1.361 second
Computing log-likelihood of 98 initial trees ... 4.530 seconds
Current best score: -9885.203

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -9871.310
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 4: -9867.255
Iteration 10 / LogL: -9879.840 / Time: 0h:0m:16s
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 12: -9865.050
Iteration 20 / LogL: -9890.182 / Time: 0h:0m:22s
Finish initializing candidate tree set (20)
Current best tree score: -9865.050 / CPU time: 18.763
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 22: -9863.685
Iteration 30 / LogL: -9873.571 / Time: 0h:0m:31s (0h:1m:39s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 35: -9863.532
Iteration 40 / LogL: -9916.789 / Time: 0h:0m:39s (0h:1m:36s left)
Iteration 50 / LogL: -9865.110 / Time: 0h:0m:47s (0h:1m:22s left)
Iteration 60 / LogL: -9925.588 / Time: 0h:0m:55s (0h:1m:10s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 64: -9863.422
Iteration 70 / LogL: -9894.578 / Time: 0h:1m:4s (0h:1m:27s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 74: -9863.389
Iteration 80 / LogL: -9896.740 / Time: 0h:1m:13s (0h:1m:26s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 85: -9862.640
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 88: -9861.140
Iteration 90 / LogL: -9863.487 / Time: 0h:1m:22s (0h:1m:30s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 100: -9860.371
Iteration 100 / LogL: -9860.371 / Time: 0h:1m:29s (0h:1m:30s left)
Iteration 110 / LogL: -9861.132 / Time: 0h:1m:38s (0h:1m:21s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 113: -9860.055
Iteration 120 / LogL: -9868.215 / Time: 0h:1m:45s (0h:1m:22s left)
Iteration 130 / LogL: -9860.568 / Time: 0h:1m:53s (0h:1m:13s left)
BETTER TREE FOUND at iteration 136: -9860.053
BETTER TREE FOUND at iteration 140: -9860.046
Iteration 140 / LogL: -9860.046 / Time: 0h:2m:1s (0h:1m:27s left)
Iteration 150 / LogL: -9861.621 / Time: 0h:2m:8s (0h:1m:17s left)
Iteration 160 / LogL: -9861.759 / Time: 0h:2m:16s (0h:1m:8s left)
Iteration 170 / LogL: -9877.400 / Time: 0h:2m:24s (0h:0m:59s left)
BETTER TREE FOUND at iteration 177: -9860.046
Iteration 180 / LogL: -9860.059 / Time: 0h:2m:33s (0h:1m:23s left)
Iteration 190 / LogL: -9861.036 / Time: 0h:2m:41s (0h:1m:14s left)
Iteration 200 / LogL: -9877.283 / Time: 0h:2m:49s (0h:1m:5s left)
Iteration 210 / LogL: -9875.817 / Time: 0h:2m:56s (0h:0m:56s left)
Iteration 220 / LogL: -9860.854 / Time: 0h:3m:4s (0h:0m:48s left)
Iteration 230 / LogL: -9860.059 / Time: 0h:3m:12s (0h:0m:39s left)
Iteration 240 / LogL: -9890.753 / Time: 0h:3m:20s (0h:0m:31s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 244: -9859.971
Iteration 250 / LogL: -9860.055 / Time: 0h:3m:29s (0h:1m:18s left)
Iteration 260 / LogL: -9860.335 / Time: 0h:3m:36s (0h:1m:10s left)
Iteration 270 / LogL: -9892.312 / Time: 0h:3m:44s (0h:1m:1s left)
Iteration 280 / LogL: -9868.546 / Time: 0h:3m:52s (0h:0m:53s left)
Iteration 290 / LogL: -9860.055 / Time: 0h:4m:0s (0h:0m:44s left)
Iteration 300 / LogL: -9868.638 / Time: 0h:4m:8s (0h:0m:36s left)
Iteration 310 / LogL: -9886.717 / Time: 0h:4m:16s (0h:0m:28s left)
Iteration 320 / LogL: -9869.503 / Time: 0h:4m:23s (0h:0m:19s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 329: -9859.786
Iteration 330 / LogL: -9861.327 / Time: 0h:4m:31s (0h:1m:21s left)
BETTER TREE FOUND at iteration 336: -9859.784
Iteration 340 / LogL: -9859.968 / Time: 0h:4m:39s (0h:1m:19s left)
BETTER TREE FOUND at iteration 345: -9859.782
Iteration 350 / LogL: -9955.081 / Time: 0h:4m:46s (0h:1m:17s left)
Iteration 360 / LogL: -9897.705 / Time: 0h:4m:53s (0h:1m:9s left)
Iteration 370 / LogL: -9880.887 / Time: 0h:5m:1s (0h:1m:1s left)
Iteration 380 / LogL: -9870.859 / Time: 0h:5m:9s (0h:0m:53s left)
BETTER TREE FOUND at iteration 390: -9859.781
Iteration 390 / LogL: -9859.781 / Time: 0h:5m:17s (0h:1m:21s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 398: -9859.743
Iteration 400 / LogL: -9893.460 / Time: 0h:5m:25s (0h:1m:19s left)
Iteration 410 / LogL: -9860.116 / Time: 0h:5m:32s (0h:1m:11s left)
Iteration 420 / LogL: -9859.789 / Time: 0h:5m:40s (0h:1m:3s left)
Iteration 430 / LogL: -9878.423 / Time: 0h:5m:48s (0h:0m:55s left)
Iteration 440 / LogL: -9861.385 / Time: 0h:5m:56s (0h:0m:47s left)
BETTER TREE FOUND at iteration 445: -9859.738
Iteration 450 / LogL: -9940.176 / Time: 0h:6m:3s (0h:1m:16s left)
Iteration 460 / LogL: -9872.571 / Time: 0h:6m:11s (0h:1m:8s left)
Iteration 470 / LogL: -9867.867 / Time: 0h:6m:19s (0h:1m:0s left)
Iteration 480 / LogL: -9859.786 / Time: 0h:6m:27s (0h:0m:52s left)
Iteration 490 / LogL: -9859.738 / Time: 0h:6m:35s (0h:0m:44s left)
Iteration 500 / LogL: -9860.038 / Time: 0h:6m:42s (0h:0m:36s left)
Iteration 510 / LogL: -9890.951 / Time: 0h:6m:50s (0h:0m:28s left)
Iteration 520 / LogL: -9859.926 / Time: 0h:6m:57s (0h:0m:20s left)
Iteration 530 / LogL: -9861.604 / Time: 0h:7m:4s (0h:0m:12s left)
Iteration 540 / LogL: -9861.051 / Time: 0h:7m:13s (0h:0m:4s left)
TREE SEARCH COMPLETED AFTER 546 ITERATIONS / Time: 0h:7m:17s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -9859.738
Optimal log-likelihood: -9859.736
Rate parameters:  A-C: 1.00000  A-G: 2.79051  A-T: 1.00000  C-G: 1.00000  C-T: 5.33648  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Site proportion and rates:  (0.838,0.149) (0.127,2.903) (0.035,14.416)
Parameters optimization took 1 rounds (0.141 sec)
BEST SCORE FOUND : -9859.736
Total tree length: 0.892

Total number of iterations: 546
CPU time used for tree search: 433.239 sec (0h:7m:13s)
Wall-clock time used for tree search: 433.928 sec (0h:7m:13s)
Total CPU time used: 437.219 sec (0h:7m:17s)
Total wall-clock time used: 437.928 sec (0h:7m:17s)

Analysis results written to: 
  IQ-TREE report:                18S-aligned-mafft.fasta.iqtree
  Maximum-likelihood tree:       18S-aligned-mafft.fasta.treefile
  Likelihood distances:          18S-aligned-mafft.fasta.mldist
  Screen log file:               18S-aligned-mafft.fasta.log
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/18S-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies/Gene_Trees

#### 28S Gene Tree
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s 28S-aligned-mafft.fasta>
-running 28S Mafft alignment
Output:
Reading alignment file 28S-aligned-mafft.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
Alignment has 169 sequences with 2441 columns, 896 distinct patterns
498 parsimony-informative, 220 singleton sites, 1723 constant sites
                          Gap/Ambiguity  Composition  p-value
   1  Zalmoxida                  16.96%    passed     48.55%
   2  Agoristenidae_DNA105839    30.36%    passed     99.18%
   3  Stygnopsidae_DNA104855      3.48%    passed     83.68%
   4  Assamiidae_DNA104857       52.48%    passed     11.18%
   5  Kimula                     80.25%    passed      8.01%
   6  Paktongius                  1.72%    passed     95.68%
   7  Assamiidae_DNA104859       20.16%    passed     80.02%
   8  Guasinia                   52.52%    passed      6.60%
   9  Icaleptes_DNA104053         1.47%    passed     96.05%
  10  Icaleptes                   1.52%    passed     97.06%
  11  Fissiphallius_DNA104055     1.47%    passed     88.71%
  12  Fissiphallius               1.47%    passed     89.90%
  13  Zalmoxis                    1.56%    passed     92.60%
  14  Ethobunus                   2.74%    passed     87.46%
  15  Stenostygnus_DNA104847      1.47%    passed     90.13%
  16  Stenostygnus_DNA104850     19.01%    passed     78.16%
  17  Stenostygnus_DNA104848     20.16%    passed     83.59%
  18  Stenostygnus_DNA104849     20.16%    passed     83.77%
  19  Fijicolana                 19.50%    passed     71.49%
  20  Pellobunus                  1.64%    passed     98.91%
  21  Glysterus                   1.52%    passed     94.40%
  22  Heterocranaus               1.52%    passed     95.02%
  23  Goniosoma                   3.32%    passed     96.49%
  24  Megapachylus                1.88%    passed     97.16%
  25  Zygopachylus                1.60%    passed     96.94%
  26  Stygnoplus                  1.52%    passed     93.16%
  27  IC_DNA103729                8.32%    passed     65.76%
  28  IC_DNA104071                1.80%    passed     90.84%
  29  IC_DNA104070                1.47%    passed     95.52%
  30  Karos                       1.47%    passed     98.37%
  31  Metalibitia                 1.64%    passed     93.81%
  32  Icaleptes_DNA104056         3.32%    passed     87.45%
  33  Metabiantes_DNA100703       1.47%    passed     98.58%
  34  Metabiantes_DNA100704       4.14%    passed     93.74%
  35  Cynortula                   1.97%    passed     96.12%
  36  Dongmoa                     1.47%    passed     82.07%
  37  Lomanius_DNA104935          1.52%    passed     71.95%
  38  Santobius_DNA104930         3.20%    passed     91.93%
  39  Santobius_DNA104931         3.20%    passed     93.28%
  40  Hoplobunus                 17.94%    passed     62.50%
  41  Stygnopsidae_DNA104856      1.47%    passed     94.87%
  42  Caenoncopus                 1.47%    passed     61.51%
  43  Sandokan                    1.47%    passed     68.52%
  44  Gnomulus                    1.52%    passed     89.51%
  45  Palaeoncopus                1.60%    passed     76.55%
  46  Lomanius_DNA104934          3.93%    passed     58.67%
  47  Conomma                     1.47%    passed     96.09%
  48  Mitraceras                  1.52%    passed     98.42%
  49  Haasus_sp                  22.41%    passed     76.90%
  50  Haasus_judaeus             49.04%    passed     15.67%
  51  IC_DNA102668                1.47%    passed     66.22%
  52  IC_DNA102669                3.24%    passed     71.04%
  53  Fissiphallius_DNA104057     3.36%    passed     80.98%
  54  Stygnommatidae_DNA105636    3.36%    passed     87.63%
  55  Epedanidae_DNA104861        1.56%    passed     98.92%
  56  Epedanidae_DNA104862        1.76%    passed     99.30%
  57  Trionyxella                 1.47%    passed     98.49%
  58  Equitius                    1.60%    passed     83.36%
  59  Triaenobunus                1.47%    passed     87.89%
  60  Larifuga                    1.52%    passed     81.38%
  61  Rostromontia               49.00%    passed     50.61%
  62  Peltonychia                 1.47%    passed     96.26%
  63  Santinezia                  5.00%    passed     92.35%
  64  Escadabius                  3.40%    passed     83.00%
  65  Arulla_DNA102666           22.98%    passed     92.74%
  66  Assamiidae_DNA104069       21.38%    passed     90.40%
  67  Paraselenca                21.22%    passed     65.88%
  68  Jarmilana                  67.60%    passed     28.76%
  69  Minuella                   12.33%    passed     74.38%
  70  Erebomaster                49.61%    passed     35.38%
  71  IC_DNA103572                9.05%    passed     39.73%
  72  As009_Cam                  66.49%    passed     40.99%
  73  As045_Cam                  66.57%    passed     41.61%
  74  As011_Gab                  66.49%    passed     44.35%
  75  As025_Gab                  66.90%    passed     49.64%
  76  As056_Gab                  66.53%    passed     50.43%
  77  As069_Lib                  66.65%    passed     49.65%
  78  As027_Gab                  66.49%    passed     40.38%
  79  As018_Gab                  66.49%    passed     44.35%
  80  As021_Gab                  66.49%    passed     51.71%
  81  As104_EAus                 66.53%    passed     41.28%
  82  As032_Gab                  66.53%    passed     46.83%
  83  As039_Gab                  66.49%    passed     44.64%
  84  As059_Gab                  66.49%    passed     37.00%
  85  As070_Lib                  67.14%    passed     38.90%
  86  As058_Gab                  66.53%    passed     38.13%
  87  As101_NAus                 66.57%    passed     47.28%
  88  As040_Gab                  66.49%    passed     36.40%
  89  As031_Gab                  66.49%    passed     42.36%
  90  As103_NAus                 67.68%    passed     25.96%
  91  As105_Phil                 66.49%    passed     26.37%
  92  As028_Gab                  66.90%    passed     53.70%
  93  As041_Gab                  66.49%    passed     52.29%
  94  As014_Gab                  68.46%    passed     27.48%
  95  As080_Indo                 66.57%    passed     69.39%
  96  As081_Indo                 66.78%    passed     66.10%
  97  As082_Cam                  66.65%    passed     69.51%
  98  As090_EAus                 66.49%    passed     72.59%
  99  As087_WAus                 66.53%    passed     76.97%
 100  As096_PNG                  66.49%    passed     79.25%
 101  As097_PNG                  66.73%    passed     78.91%
 102  As094_PNG                  66.65%    passed     70.19%
 103  As098_Phil                 66.53%    passed     71.45%
 104  As109_Thai                 66.82%    passed     56.37%
 105  As121_Laos                 66.94%    passed     54.17%
 106  As110_Laos                 66.90%    passed     57.58%
 107  As120_Thai                 66.49%    passed     52.29%
 108  As132_Mala                 66.57%    passed     57.90%
 109  As085_WAus                 66.57%    passed     69.71%
 110  As083_Cam                  68.21%    passed     46.08%
 111  As127_Laos                 67.31%    passed     54.01%
 112  As134_Viet                 66.49%    passed     66.02%
 113  As133_Indo                 66.86%    passed     69.25%
 114  As117_Thai                 66.53%    passed     77.67%
 115  As020_Gab                  66.49%    passed     24.05%
 116  As026_Gab                  66.49%    passed     21.76%
 117  As038_Gab                  66.53%    passed     23.06%
 118  As135_Viet                 66.98%    passed     40.14%
 119  Baculigerus                 8.56%    passed     97.50%
 120  Pseudoepedanus              8.60%    passed     99.49%
 121  Stygnomma                  26.59%    passed     62.84%
 122  Biantidae_DNA105668        27.73%    passed     74.35%
 123  Triaenonychidae             6.84%    passed     98.92%
 124  As092_EAus                 66.57%    passed     32.86%
 125  As099_Dibunus              65.22%    passed     45.15%
 126  Op104_Dibuninae            66.24%    passed     35.81%
 127  Op105_Dibunus              66.37%    passed     22.36%
 128  Op106_Toccolus             66.12%    passed     44.60%
 129  Op107_Nanepedanus_rufus    66.20%    passed     35.82%
 130  Op049_Beloniscus           66.33%    passed     69.35%
 131  Op103_Bupares              66.49%    passed     44.38%
 132  Stygnopsidae_DNA103882     29.25%    passed     97.13%
 133  Samoidae                   29.54%    passed     97.10%
 134  Theromaster                 5.57%    passed     96.63%
 135  Bunofagea                   4.63%    passed     87.00%
 136  Scotolemon                  1.52%    passed     92.46%
 137  Scotolemon_lespesi         49.04%    passed     73.02%
 138  Holoscotolemon             24.91%    passed     70.05%
 139  Trojanella                 56.12%    passed     13.32%
 140  As106_Thai                 66.94%    passed     72.58%
 141  As129_Laos                 80.34%    passed     63.80%
 142  As115_Thai                 66.57%    passed     67.38%
 143  As126_Laos                 67.02%    passed     65.96%
 144  As122_Laos                 67.35%    passed     58.54%
 145  As140_Viet                 67.23%    passed     57.10%
 146  As136_Viet                 66.61%    passed     66.39%
 147  As123_Viet                 67.23%    passed     60.03%
 148  As118_Thai                 66.69%    passed     73.99%
 149  As131_Mala                 66.78%    passed     59.09%
 150  Epedanidae_DNA104062       49.04%    passed     29.14%
 151  Epedanidae_DNA104066       29.29%    passed     93.72%
 152  Tithaeus                   29.86%    passed     98.90%
 153  Epedanidae_DNA104068       39.90%    passed     84.38%
 154  Chilon                     67.39%    passed     67.39%
 155  Metabiantes_DNA100335      34.08%    passed     26.77%
 156  Maiorerus                  87.14%    passed     45.32%
 157  Zuma                       88.00%    passed     46.98%
 158  Martensiellus               6.35%    passed     76.07%
 159  Synthetonychia             87.96%    passed     74.83%
 160  Bishopella                 34.62%    passed     97.64%
 161  Vonones                    87.18%    passed     16.55%
 162  Caddo                      87.14%    passed     36.22%
 163  Pantopsalis                49.12%    passed     70.86%
 164  Protolophus                 1.64%    passed     11.88%
 165  Dendrolasma                 2.13%    failed      0.49%
 166  Trogulus                    8.93%    failed      0.02%
 167  Hesperonemastoma            4.88%    failed      0.07%
 168  Troglosiro                 31.63%    passed     92.46%
 169  Fumontana                  40.39%    passed     78.98%
WARNING: 76 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      37.15%  3 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.014 seconds
NOTE: ModelFinder requires 49 MB RAM!
ModelFinder will test 286 DNA models (sample size: 2441) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  JC            18910.325    335 38490.649    38597.595    40433.704
  2  JC+I          17417.463    336 35506.926    35614.561    37455.781
  3  JC+G4         17005.312    336 34682.625    34790.260    36631.480
  4  JC+I+G4       16819.222    337 34312.444    34420.771    36267.099
  5  JC+R2         16929.835    337 34533.670    34641.998    36488.325
  6  JC+R3         16720.199    339 34118.397    34228.117    36084.653
  7  JC+R4         16690.355    341 34062.710    34173.831    36040.565
  8  JC+R5         16678.827    343 34043.653    34156.187    36033.109
  9  JC+R6         16667.728    345 34025.455    34139.412    36026.511
 10  JC+R7         16664.218    347 34022.437    34137.827    36035.093
 14  F81+F         18963.699    338 38603.399    38712.421    40563.854
 15  F81+F+I       17464.388    339 35606.777    35716.496    37573.032
 16  F81+F+G4      17059.543    339 34797.085    34906.804    36763.340
 17  F81+F+I+G4    16877.293    340 34434.586    34545.005    36406.641
 18  F81+F+R2      16984.521    340 34649.042    34759.461    36621.098
 19  F81+F+R3      16779.599    342 34243.198    34355.025    36226.854
 20  F81+F+R4      16750.385    344 34188.770    34302.014    36184.026
 21  F81+F+R5      16739.164    346 34170.327    34285.000    36177.184
 22  F81+F+R6      16728.416    348 34152.831    34268.942    36171.288
 23  F81+F+R7      16725.059    350 34150.119    34267.679    36180.176
 27  K2P           18614.088    336 37900.177    38007.812    39849.032
 28  K2P+I         17109.109    337 34892.218    35000.545    36846.873
 29  K2P+G4        16685.788    337 34045.575    34153.903    36000.230
 30  K2P+I+G4      16494.188    338 33664.376    33773.398    35624.832
 31  K2P+R2        16611.375    338 33898.751    34007.773    35859.206
 32  K2P+R3        16395.959    340 33471.919    33582.338    35443.974
 33  K2P+R4        16364.533    342 33413.065    33524.892    35396.721
 34  K2P+R5        16351.783    344 33391.565    33504.809    35386.821
 35  K2P+R6        16339.270    346 33370.539    33485.212    35377.396
 36  K2P+R7        16335.211    348 33366.423    33482.534    35384.880
 40  HKY+F         18663.674    339 38005.347    38115.067    39971.603
 41  HKY+F+I       17131.445    340 34942.890    35053.309    36914.946
 42  HKY+F+G4      16711.259    340 34102.517    34212.936    36074.573
 43  HKY+F+I+G4    16522.185    341 33726.369    33837.491    35704.225
 44  HKY+F+R2      16639.883    341 33961.766    34072.887    35939.622
 45  HKY+F+R3      16430.609    343 33547.217    33659.751    35536.673
 46  HKY+F+R4      16399.153    345 33488.306    33602.263    35489.362
 47  HKY+F+R5      16386.210    347 33466.421    33581.811    35479.077
 48  HKY+F+R6      16373.388    349 33444.776    33561.610    35469.033
 49  HKY+F+R7      16369.187    351 33440.374    33558.662    35476.231
 53  TNe           18543.377    337 37760.754    37869.082    39715.409
 54  TNe+I         17076.502    338 34829.003    34938.025    36789.458
 55  TNe+G4        16659.045    338 33994.090    34103.112    35954.545
 56  TNe+I+G4      16471.581    339 33621.161    33730.880    35587.416
 57  TNe+R2        16584.534    339 33847.068    33956.787    35813.324
 58  TNe+R3        16369.352    341 33420.703    33531.825    35398.559
 59  TNe+R4        16337.947    343 33361.894    33474.429    35351.350
 60  TNe+R5        16324.936    345 33339.872    33453.829    35340.928
 61  TNe+R6        16311.964    347 33317.928    33433.319    35330.585
 62  TNe+R7        16307.724    349 33313.448    33430.282    35337.705
 66  TN+F          18552.656    340 37785.311    37895.730    39757.367
 67  TN+F+I        17075.344    341 34832.688    34943.809    36810.543
 68  TN+F+G4       16651.578    341 33985.157    34096.278    35963.012
 69  TN+F+I+G4     16472.201    342 33628.403    33740.229    35612.059
 70  TN+F+R2       16583.394    342 33850.787    33962.614    35834.443
 71  TN+F+R3       16373.763    344 33435.526    33548.770    35430.782
 72  TN+F+R4       16342.195    346 33376.390    33491.062    35383.246
 73  TN+F+R5       16328.640    348 33353.280    33469.391    35371.737
 74  TN+F+R6       16314.940    350 33329.881    33447.441    35359.938
 75  TN+F+R7       16310.396    352 33324.793    33443.812    35366.450
 79  K3P           18613.869    337 37901.738    38010.065    39856.393
 80  K3P+I         17108.869    338 34893.738    35002.760    36854.193
 81  K3P+G4        16685.564    338 34047.128    34156.150    36007.584
 82  K3P+I+G4      16493.854    339 33665.709    33775.428    35631.964
 83  K3P+R2        16611.145    339 33900.289    34010.008    35866.544
 84  K3P+R3        16395.525    341 33473.051    33584.172    35450.906
 85  K3P+R4        16363.974    343 33413.948    33526.482    35403.404
 86  K3P+R5        16351.116    345 33392.232    33506.189    35393.288
 87  K3P+R6        16338.490    347 33370.981    33486.371    35383.637
 88  K3P+R7        16334.398    349 33366.795    33483.629    35391.052
 92  K3Pu+F        18662.990    340 38005.980    38116.399    39978.035
 93  K3Pu+F+I      17130.118    341 34942.236    35053.357    36920.091
 94  K3Pu+F+G4     16709.794    341 34101.589    34212.710    36079.445
 95  K3Pu+F+I+G4   16520.587    342 33725.174    33837.000    35708.829
 96  K3Pu+F+R2     16638.506    342 33961.013    34072.839    35944.668
 97  K3Pu+F+R3     16429.042    344 33546.085    33659.329    35541.341
 98  K3Pu+F+R4     16397.476    346 33486.953    33601.625    35493.809
 99  K3Pu+F+R5     16384.381    348 33464.761    33580.872    35483.218
100  K3Pu+F+R6     16371.407    350 33442.814    33560.374    35472.871
101  K3Pu+F+R7     16367.174    352 33438.348    33557.367    35480.006
105  TPM2+F        18660.246    340 38000.493    38110.912    39972.548
106  TPM2+F+I      17117.788    341 34917.575    35028.697    36895.431
107  TPM2+F+G4     16700.169    341 34082.338    34193.459    36060.193
108  TPM2+F+I+G4   16512.178    342 33708.357    33820.183    35692.013
109  TPM2+F+R2     16631.236    342 33946.472    34058.299    35930.128
110  TPM2+F+R3     16423.452    344 33534.904    33648.148    35530.160
111  TPM2+F+R4     16391.981    346 33475.961    33590.634    35482.818
112  TPM2+F+R5     16378.836    348 33453.672    33569.783    35472.128
113  TPM2+F+R6     16365.792    350 33431.583    33549.143    35461.640
114  TPM2+F+R7     16361.523    352 33427.046    33546.065    35468.703
118  TPM2u+F       18660.252    340 38000.504    38110.923    39972.559
119  TPM2u+F+I     17117.796    341 34917.592    35028.713    36895.447
120  TPM2u+F+G4    16700.172    341 34082.344    34193.465    36060.199
121  TPM2u+F+I+G4  16512.191    342 33708.381    33820.208    35692.037
122  TPM2u+F+R2    16631.259    342 33946.518    34058.345    35930.174
123  TPM2u+F+R3    16423.458    344 33534.915    33648.160    35530.171
124  TPM2u+F+R4    16391.983    346 33475.966    33590.639    35482.823
125  TPM2u+F+R5    16378.829    348 33453.658    33569.769    35472.115
126  TPM2u+F+R6    16365.777    350 33431.554    33549.114    35461.611
127  TPM2u+F+R7    16361.507    352 33427.014    33546.033    35468.671
131  TPM3+F        18656.213    340 37992.427    38102.846    39964.482
132  TPM3+F+I      17120.116    341 34922.231    35033.353    36900.087
133  TPM3+F+G4     16691.439    341 34064.878    34175.999    36042.733
134  TPM3+F+I+G4   16504.636    342 33693.271    33805.098    35676.927
135  TPM3+F+R2     16620.773    342 33925.546    34037.373    35909.202
136  TPM3+F+R3     16414.363    344 33516.726    33629.970    35511.982
137  TPM3+F+R4     16382.805    346 33457.611    33572.283    35464.467
138  TPM3+F+R5     16369.686    348 33435.372    33551.483    35453.829
139  TPM3+F+R6     16356.634    350 33413.267    33530.827    35443.325
140  TPM3+F+R7     16352.359    352 33408.718    33527.737    35450.376
144  TPM3u+F       18656.216    340 37992.432    38102.851    39964.487
145  TPM3u+F+I     17120.117    341 34922.234    35033.355    36900.090
146  TPM3u+F+G4    16691.440    341 34064.880    34176.002    36042.736
147  TPM3u+F+I+G4  16504.664    342 33693.329    33805.155    35676.984
148  TPM3u+F+R2    16620.764    342 33925.529    34037.355    35909.185
149  TPM3u+F+R3    16414.364    344 33516.727    33629.971    35511.983
150  TPM3u+F+R4    16382.806    346 33457.612    33572.284    35464.468
151  TPM3u+F+R5    16369.681    348 33435.363    33551.474    35453.820
152  TPM3u+F+R6    16356.628    350 33413.256    33530.816    35443.313
153  TPM3u+F+R7    16352.357    352 33408.714    33527.733    35450.371
157  TIMe          18543.187    338 37762.374    37871.396    39722.829
158  TIMe+I        17076.289    339 34830.578    34940.297    36796.833
159  TIMe+G4       16658.805    339 33995.611    34105.330    35961.866
160  TIMe+I+G4     16471.256    340 33622.512    33732.931    35594.568
161  TIMe+R2       16584.299    340 33848.598    33959.018    35820.654
162  TIMe+R3       16368.960    342 33421.921    33533.747    35405.577
163  TIMe+R4       16337.404    344 33362.807    33476.051    35358.063
164  TIMe+R5       16324.333    346 33340.667    33455.339    35347.523
165  TIMe+R6       16311.333    348 33318.665    33434.776    35337.122
166  TIMe+R7       16307.082    350 33314.164    33431.724    35344.221
170  TIM+F         18551.955    341 37785.911    37897.032    39763.766
171  TIM+F+I       17074.019    342 34832.038    34943.865    36815.694
172  TIM+F+G4      16649.956    342 33983.912    34095.739    35967.568
173  TIM+F+I+G4    16470.637    343 33627.273    33739.807    35616.729
174  TIM+F+R2      16581.926    343 33849.852    33962.386    35839.308
175  TIM+F+R3      16372.198    345 33434.395    33548.352    35435.451
176  TIM+F+R4      16340.444    347 33374.889    33490.279    35387.545
177  TIM+F+R5      16326.840    349 33351.681    33468.515    35375.937
178  TIM+F+R6      16313.062    351 33328.124    33446.413    35363.982
179  TIM+F+R7      16308.491    353 33322.981    33442.734    35370.439
183  TIM2e         18543.361    338 37762.722    37871.743    39723.177
184  TIM2e+I       17073.388    339 34824.777    34934.496    36791.032
185  TIM2e+G4      16657.693    339 33993.385    34103.104    35959.641
186  TIM2e+I+G4    16470.980    340 33621.960    33732.379    35594.016
187  TIM2e+R2      16583.969    340 33847.939    33958.358    35819.994
188  TIM2e+R3      16369.167    342 33422.333    33534.160    35405.989
189  TIM2e+R4      16337.653    344 33363.306    33476.551    35358.562
190  TIM2e+R5      16324.602    346 33341.203    33455.876    35348.060
191  TIM2e+R6      16311.580    348 33319.159    33435.270    35337.616
192  TIM2e+R7      16307.318    350 33314.635    33432.195    35344.693
196  TIM2+F        18549.275    341 37780.549    37891.671    39758.405
197  TIM2+F+I      17061.657    342 34807.313    34919.140    36790.969
198  TIM2+F+G4     16636.030    342 33956.060    34067.887    35939.716
199  TIM2+F+I+G4   16460.631    343 33607.262    33719.796    35596.718
200  TIM2+F+R2     16572.713    343 33831.427    33943.961    35820.883
201  TIM2+F+R3     16365.382    345 33420.764    33534.721    35421.820
202  TIM2+F+R4     16333.995    347 33361.990    33477.380    35374.646
203  TIM2+F+R5     16320.779    349 33339.558    33456.393    35363.815
204  TIM2+F+R6     16306.690    351 33315.379    33433.667    35351.236
205  TIM2+F+R7     16302.080    353 33310.160    33429.913    35357.617
209  TIM3e         18541.105    338 37758.209    37867.231    39718.664
210  TIM3e+I       17071.895    339 34821.789    34931.508    36788.045
211  TIM3e+G4      16652.638    339 33983.276    34092.995    35949.531
212  TIM3e+I+G4    16466.115    340 33612.230    33722.649    35584.286
213  TIM3e+R2      16577.723    340 33835.446    33945.865    35807.502
214  TIM3e+R3      16363.962    342 33411.924    33523.750    35395.579
215  TIM3e+R4      16332.278    344 33352.555    33465.800    35347.812
216  TIM3e+R5      16319.179    346 33330.358    33445.030    35337.214
217  TIM3e+R6      16306.156    348 33308.312    33424.423    35326.768
218  TIM3e+R7      16301.915    350 33303.830    33421.389    35333.887
222  TIM3+F        18545.078    341 37772.157    37883.278    39750.012
223  TIM3+F+I      17064.169    342 34812.338    34924.164    36795.994
224  TIM3+F+G4     16634.033    342 33952.067    34063.893    35935.722
225  TIM3+F+I+G4   16456.532    343 33599.065    33711.599    35588.521
226  TIM3+F+R2     16566.724    343 33819.448    33931.982    35808.904
227  TIM3+F+R3     16358.661    345 33407.321    33521.278    35408.378
228  TIM3+F+R4     16326.518    347 33347.035    33462.426    35359.692
229  TIM3+F+R5     16313.121    349 33324.242    33441.076    35348.499
230  TIM3+F+R6     16299.237    351 33300.474    33418.762    35336.331
231  TIM3+F+R7     16294.754    353 33295.507    33415.260    35342.965
235  TVMe          18611.706    339 37901.412    38011.131    39867.667
236  TVMe+I        17101.113    340 34882.226    34992.645    36854.281
237  TVMe+G4       16677.500    340 34035.000    34145.419    36007.055
238  TVMe+I+G4     16487.082    341 33656.163    33767.285    35634.019
239  TVMe+R2       16603.073    341 33888.146    33999.268    35866.002
240  TVMe+R3       16389.576    343 33465.152    33577.687    35454.608
241  TVMe+R4       16357.978    345 33405.955    33519.912    35407.011
242  TVMe+R5       16345.306    347 33384.611    33500.002    35397.268
243  TVMe+R6       16332.299    349 33362.599    33479.433    35386.856
244  TVMe+R7       16328.215    351 33358.430    33476.718    35394.287
248  TVM+F         18652.551    342 37989.103    38100.929    39972.758
249  TVM+F+I       17106.210    343 34898.420    35010.954    36887.876
250  TVM+F+G4      16680.076    343 34046.151    34158.685    36035.607
251  TVM+F+I+G4    16494.378    344 33676.756    33790.000    35672.012
252  TVM+F+R2      16611.802    344 33911.604    34024.848    35906.860
253  TVM+F+R3      16406.867    346 33505.734    33620.406    35512.590
254  TVM+F+R4      16375.358    348 33446.716    33562.827    35465.173
255  TVM+F+R5      16362.576    350 33425.152    33542.711    35455.209
256  TVM+F+R6      16348.849    352 33401.698    33520.717    35443.356
257  TVM+F+R7      16344.471    354 33396.942    33517.431    35450.200
261  SYM           18540.984    340 37761.969    37872.388    39734.024
262  SYM+I         17068.687    341 34819.374    34930.495    36797.229
263  SYM+G4        16651.225    341 33984.450    34095.572    35962.306
264  SYM+I+G4      16465.497    342 33614.994    33726.820    35598.650
265  SYM+R2        16577.086    342 33838.172    33949.998    35821.827
266  SYM+R3        16363.630    344 33415.259    33528.504    35410.515
267  SYM+R4        16331.945    346 33355.890    33470.562    35362.746
268  SYM+R5        16318.828    348 33333.655    33449.766    35352.112
269  SYM+R6        16305.781    350 33311.561    33429.121    35341.618
270  SYM+R7        16301.512    352 33307.023    33426.042    35348.681
274  GTR+F         18541.495    343 37768.991    37881.525    39758.447
275  GTR+F+I       17050.270    344 34788.541    34901.785    36783.797
276  GTR+F+G4      16618.178    344 33924.356    34037.600    35919.612
277  GTR+F+I+G4    16444.755    345 33579.511    33693.468    35580.567
278  GTR+F+R2      16555.740    345 33801.480    33915.437    35802.536
279  GTR+F+R3      16349.905    347 33393.809    33509.199    35406.466
280  GTR+F+R4      16318.020    349 33334.040    33450.874    35358.297
281  GTR+F+R5      16304.807    351 33311.614    33429.902    35347.471
282  GTR+F+R6      16290.649    353 33287.299    33407.051    35334.756
283  GTR+F+R7      16286.086    355 33282.172    33403.399    35341.229
Akaike Information Criterion:           GTR+F+R7
Corrected Akaike Information Criterion: GTR+F+R7
Bayesian Information Criterion:         TIM3e+R6
Best-fit model: TIM3e+R6 chosen according to BIC

All model information printed to 28S-aligned-mafft.fasta.model.gz
CPU time for ModelFinder: 180.518 seconds (0h:3m:0s)
Wall-clock time for ModelFinder: 180.991 seconds (0h:3m:0s)

NOTE: 29 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -17618.841
2. Current log-likelihood: -16454.228
3. Current log-likelihood: -16344.075
4. Current log-likelihood: -16330.040
5. Current log-likelihood: -16322.778
6. Current log-likelihood: -16317.962
7. Current log-likelihood: -16314.925
8. Current log-likelihood: -16313.089
9. Current log-likelihood: -16311.926
10. Current log-likelihood: -16311.060
11. Current log-likelihood: -16310.399
12. Current log-likelihood: -16309.879
13. Current log-likelihood: -16309.458
14. Current log-likelihood: -16309.116
15. Current log-likelihood: -16308.830
16. Current log-likelihood: -16308.593
17. Current log-likelihood: -16308.392
18. Current log-likelihood: -16308.222
19. Current log-likelihood: -16308.076
20. Current log-likelihood: -16307.951
21. Current log-likelihood: -16307.844
Optimal log-likelihood: -16307.750
Rate parameters:  A-C: 0.76603  A-G: 2.05984  A-T: 1.00000  C-G: 0.76603  C-T: 3.54669  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Site proportion and rates:  (0.803,0.147) (0.012,0.483) (0.055,2.078) (0.068,2.090) (0.041,5.768) (0.023,16.809)
Parameters optimization took 21 rounds (13.905 sec)
Computing ML distances based on estimated model parameters... 0.325 sec
WARNING: Some pairwise ML distances are too long (saturated)
Computing BIONJ tree...
0.010 seconds
Log-likelihood of BIONJ tree: -17062.915
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 1.453 second
Computing log-likelihood of 98 initial trees ... 7.844 seconds
Current best score: -16245.950

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -16212.887
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 2: -16207.896
Iteration 10 / LogL: -16230.165 / Time: 0h:0m:36s
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 18: -16204.370
Iteration 20 / LogL: -16227.239 / Time: 0h:0m:46s
Finish initializing candidate tree set (20)
Current best tree score: -16204.370 / CPU time: 32.388
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 21: -16200.019
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 25: -16198.287
BETTER TREE FOUND at iteration 30: -16198.286
Iteration 30 / LogL: -16198.286 / Time: 0h:1m:1s (0h:3m:32s left)
Iteration 40 / LogL: -16205.480 / Time: 0h:1m:14s (0h:2m:52s left)
Iteration 50 / LogL: -16397.914 / Time: 0h:1m:28s (0h:2m:23s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 53: -16197.457
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 55: -16197.220
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 60: -16196.091
Iteration 60 / LogL: -16196.091 / Time: 0h:1m:43s (0h:2m:55s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 68: -16195.229
Iteration 70 / LogL: -16207.698 / Time: 0h:1m:56s (0h:2m:45s left)
Iteration 80 / LogL: -16199.587 / Time: 0h:2m:8s (0h:2m:23s left)
Iteration 90 / LogL: -16200.947 / Time: 0h:2m:20s (0h:2m:2s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 94: -16194.331
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 98: -16191.048
Iteration 100 / LogL: -16195.507 / Time: 0h:2m:33s (0h:2m:32s left)
Iteration 110 / LogL: -16237.482 / Time: 0h:2m:46s (0h:2m:14s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 113: -16190.735
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 116: -16188.416
Iteration 120 / LogL: -16210.300 / Time: 0h:2m:59s (0h:2m:24s left)
Iteration 130 / LogL: -16191.028 / Time: 0h:3m:12s (0h:2m:8s left)
Iteration 140 / LogL: -16189.291 / Time: 0h:3m:25s (0h:1m:52s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 145: -16188.379
Iteration 150 / LogL: -16192.840 / Time: 0h:3m:38s (0h:2m:19s left)
BETTER TREE FOUND at iteration 151: -16188.378
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 155: -16188.197
Iteration 160 / LogL: -16193.414 / Time: 0h:3m:50s (0h:2m:17s left)
Iteration 170 / LogL: -16188.378 / Time: 0h:4m:3s (0h:2m:2s left)
BETTER TREE FOUND at iteration 173: -16188.192
Iteration 180 / LogL: -16191.812 / Time: 0h:4m:15s (0h:2m:12s left)
Iteration 190 / LogL: -16188.475 / Time: 0h:4m:27s (0h:1m:57s left)
Iteration 200 / LogL: -16215.958 / Time: 0h:4m:40s (0h:1m:42s left)
Iteration 210 / LogL: -16202.846 / Time: 0h:4m:52s (0h:1m:28s left)
Iteration 220 / LogL: -16202.717 / Time: 0h:5m:5s (0h:1m:13s left)
Iteration 230 / LogL: -16194.471 / Time: 0h:5m:18s (0h:0m:59s left)
BETTER TREE FOUND at iteration 238: -16188.192
Iteration 240 / LogL: -16197.078 / Time: 0h:5m:30s (0h:2m:15s left)
Iteration 250 / LogL: -16200.560 / Time: 0h:5m:43s (0h:2m:1s left)
Iteration 260 / LogL: -16216.326 / Time: 0h:5m:55s (0h:1m:47s left)
BETTER TREE FOUND at iteration 266: -16188.190
Iteration 270 / LogL: -16194.131 / Time: 0h:6m:8s (0h:2m:11s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 278: -16187.922
Iteration 280 / LogL: -16220.913 / Time: 0h:6m:22s (0h:2m:14s left)
BETTER TREE FOUND at iteration 288: -16187.919
Iteration 290 / LogL: -16191.627 / Time: 0h:6m:35s (0h:2m:13s left)
Iteration 300 / LogL: -16204.014 / Time: 0h:6m:47s (0h:1m:59s left)
Iteration 310 / LogL: -16188.212 / Time: 0h:6m:59s (0h:1m:45s left)
Iteration 320 / LogL: -16188.459 / Time: 0h:7m:12s (0h:1m:32s left)
Iteration 330 / LogL: -16189.863 / Time: 0h:7m:25s (0h:1m:18s left)
BETTER TREE FOUND at iteration 335: -16187.917
Iteration 340 / LogL: -16252.366 / Time: 0h:7m:38s (0h:2m:8s left)
BETTER TREE FOUND at iteration 343: -16187.914
Iteration 350 / LogL: -16188.105 / Time: 0h:7m:51s (0h:2m:5s left)
Iteration 360 / LogL: -16188.021 / Time: 0h:8m:3s (0h:1m:51s left)
Iteration 370 / LogL: -16194.732 / Time: 0h:8m:16s (0h:1m:38s left)
Iteration 380 / LogL: -16209.522 / Time: 0h:8m:29s (0h:1m:24s left)
Iteration 390 / LogL: -16188.291 / Time: 0h:8m:42s (0h:1m:11s left)
Iteration 400 / LogL: -16198.242 / Time: 0h:8m:54s (0h:0m:57s left)
Iteration 410 / LogL: -16194.248 / Time: 0h:9m:7s (0h:0m:44s left)
BETTER TREE FOUND at iteration 413: -16187.914
Iteration 420 / LogL: -16191.939 / Time: 0h:9m:20s (0h:2m:4s left)
Iteration 430 / LogL: -16188.479 / Time: 0h:9m:33s (0h:1m:50s left)
Iteration 440 / LogL: -16189.802 / Time: 0h:9m:45s (0h:1m:37s left)
Iteration 450 / LogL: -16198.524 / Time: 0h:9m:58s (0h:1m:23s left)
Iteration 460 / LogL: -16195.189 / Time: 0h:10m:9s (0h:1m:10s left)
Iteration 470 / LogL: -16207.830 / Time: 0h:10m:22s (0h:0m:57s left)
Iteration 480 / LogL: -16188.192 / Time: 0h:10m:34s (0h:0m:43s left)
Iteration 490 / LogL: -16203.840 / Time: 0h:10m:47s (0h:0m:30s left)
Iteration 500 / LogL: -16188.768 / Time: 0h:11m:0s (0h:0m:17s left)
Iteration 510 / LogL: -16191.263 / Time: 0h:11m:13s (0h:0m:3s left)
TREE SEARCH COMPLETED AFTER 514 ITERATIONS / Time: 0h:11m:18s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -16187.914
Optimal log-likelihood: -16187.909
Rate parameters:  A-C: 0.76396  A-G: 2.11251  A-T: 1.00000  C-G: 0.76396  C-T: 3.61991  G-T: 1.00000
Base frequencies:  A: 0.250  C: 0.250  G: 0.250  T: 0.250
Site proportion and rates:  (0.808,0.151) (0.001,0.495) (0.057,1.975) (0.071,1.976) (0.041,5.955) (0.022,17.214)
Parameters optimization took 1 rounds (0.291 sec)
BEST SCORE FOUND : -16187.909
Total tree length: 1.735

Total number of iterations: 514
CPU time used for tree search: 662.973 sec (0h:11m:2s)
Wall-clock time used for tree search: 663.727 sec (0h:11m:3s)
Total CPU time used: 677.793 sec (0h:11m:17s)
Total wall-clock time used: 678.578 sec (0h:11m:18s)

Analysis results written to: 
  IQ-TREE report:                28S-aligned-mafft.fasta.iqtree
  Maximum-likelihood tree:       28S-aligned-mafft.fasta.treefile
  Likelihood distances:          28S-aligned-mafft.fasta.mldist
  Screen log file:               28S-aligned-mafft.fasta.log

Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/28S-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies/Gene_Trees

#### COI Gene Tree
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s COI-aligned-mafft.fasta>
-running COI Mafft alignment
Output:
Reading alignment file COI-aligned-mafft.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
Alignment has 120 sequences with 657 columns, 507 distinct patterns
437 parsimony-informative, 45 singleton sites, 175 constant sites
                         Gap/Ambiguity  Composition  p-value
   1  Trojanella                 1.37%    failed      0.01%
   2  Holoscotolemon            10.35%    failed      0.00%
   3  Peltonychia                1.52%    failed      0.00%
   4  Zuma                      35.31%    failed      0.09%
   5  As009_Cam                  0.46%    failed      0.00%
   6  As050_Cam                  0.46%    failed      0.00%
   7  As018_Gab                  0.46%    failed      0.00%
   8  As011_Gab                  0.46%    failed      0.00%
   9  As056_Gab                  0.61%    failed      0.00%
  10  As027_Gab                  0.46%    failed      0.00%
  11  As084_Cam                  0.46%    failed      0.00%
  12  As057_Gab                  0.00%    failed      0.00%
  13  As058_Gab                  0.15%    failed      0.00%
  14  As040_Gab                  0.15%    failed      0.00%
  15  Scotolemon                 1.37%    failed      0.00%
  16  Bishopella                 1.37%    passed      5.50%
  17  Protolophus                1.37%    failed      0.13%
  18  Pantopsalis                1.37%    failed      1.71%
  19  Hesperonemastoma           2.28%    failed      0.00%
  20  Dendrolasma                3.65%    passed     26.51%
  21  Trogulus                   1.37%    passed     44.74%
  22  Troglosiro                 1.37%    passed      8.16%
  23  Synthetonychia             1.37%    failed      0.01%
  24  Conomma                    1.83%    failed      0.00%
  25  Guasinia                  14.16%    failed      0.00%
  26  Icaleptes_DNA104056        1.37%    failed      2.50%
  27  Ethobunus                  1.37%    failed      2.07%
  28  Zalmoxis                   1.37%    failed      0.62%
  29  Kimula                     1.37%    failed      0.03%
  30  Stygnomma                  1.37%    passed      8.67%
  31  Baculigerus                7.76%    passed      6.74%
  32  Stenostygnus_DNA104848     1.37%    passed     46.34%
  33  Stenostygnus_DNA104849     1.37%    passed     25.82%
  34  Stenostygnus_DNA104847     1.37%    passed     22.59%
  35  Stenostygnus_DNA104850     1.37%    failed      3.08%
  36  Fijicolana                 8.07%    passed      7.23%
  37  Pellobunus                 1.37%    passed     13.07%
  38  Metabiantes_DNA100704      1.37%    failed      0.11%
  39  Metabiantes_DNA100335      1.37%    failed      2.02%
  40  Biantidae_DNA105668        8.07%    failed      1.08%
  41  Metabiantes_DNA100703      1.37%    failed      3.94%
  42  Trionyxella                0.00%    failed      3.28%
  43  Paktongius                 1.37%    failed      4.43%
  44  As120_Thai                 0.00%    passed     37.90%
  45  As110_Laos                 0.46%    failed      0.18%
  46  As121_Laos                 0.00%    failed      0.26%
  47  As096_PNG                  0.00%    passed     14.64%
  48  As095_PNG                  0.00%    passed     13.28%
  49  As094_PNG                  0.00%    passed      6.50%
  50  As102_NAus                 0.00%    passed     13.12%
  51  As092_EAus                 0.15%    failed      1.91%
  52  As104_EAus                 4.57%    failed      2.40%
  53  As108_Laos                 0.00%    failed      0.01%
  54  As081_Indo                 0.00%    passed      5.40%
  55  As126_Laos                 1.22%    failed      1.61%
  56  As080_Indo                 0.00%    failed      0.07%
  57  As133_Indo                 0.00%    failed      0.01%
  58  As127_Laos                 0.46%    failed      0.00%
  59  As109_Thai                 0.00%    passed     10.35%
  60  As087_WAus                 0.00%    failed      0.00%
  61  Heterocranaus              1.83%    failed      2.66%
  62  Megapachylus               2.28%    passed      5.39%
  63  Goniosoma                  1.83%    failed      1.83%
  64  Glysterus                  8.22%    passed      5.90%
  65  Agoristenidae_DNA105839    1.83%    failed      3.38%
  66  Caenoncopus                1.37%    passed     40.13%
  67  Gnomulus                   1.37%    passed     58.36%
  68  Palaeoncopus               1.37%    failed      0.15%
  69  Martensiellus              1.37%    passed     28.73%
  70  Sandokan                   1.37%    failed      3.32%
  71  Epedanidae_DNA104066       3.50%    passed     91.20%
  72  Tithaeus                   7.76%    failed      0.23%
  73  As106_Thai                 0.00%    failed      0.03%
  74  As118_Thai                 0.00%    failed      0.01%
  75  As115_Thai                 0.00%    failed      0.01%
  76  As131_Mala                 0.00%    passed     11.23%
  77  As114_Thai                 0.00%    failed      3.20%
  78  As129_Laos                 1.07%    failed      2.47%
  79  As122_Laos                 0.00%    failed      2.11%
  80  As136_Viet                 0.30%    passed     13.36%
  81  As123_Viet                 0.00%    failed      0.01%
  82  As140_Viet                 0.00%    passed     12.24%
  83  Epedanidae_DNA104062       1.37%    passed     38.87%
  84  Epedanidae_DNA104068       1.37%    passed     47.73%
  85  Zygopachylus               1.83%    passed     63.18%
  86  Cynortula                  1.98%    passed     37.28%
  87  As026_Gab                  0.91%    failed      0.00%
  88  Hoplobunus                 2.28%    passed     27.46%
  89  Stygnopsidae_DNA103882     8.68%    failed      0.00%
  90  Karos                      2.28%    failed      3.55%
  91  Stygnopsidae_DNA104855     2.28%    passed     32.69%
  92  Stygnopsidae_DNA104856     2.28%    passed     37.72%
  93  As020_Gab                  0.46%    passed      6.45%
  94  As041_Gab                  0.46%    failed      0.68%
  95  As030_Gab                  0.46%    passed      7.08%
  96  Jarmilana                  4.11%    failed      2.92%
  97  Arulla_DNA102666           0.00%    failed      0.00%
  98  As059_Gab                  1.52%    failed      0.00%
  99  IC_DNA104070               7.76%    failed      0.31%
 100  Zalmoxida                  0.00%    failed      0.01%
 101  As099_Dibunus              0.00%    failed      4.04%
 102  Op105_Dibunus              3.65%    failed      1.71%
 103  Op106_Toccolus             3.65%    failed      0.84%
 104  Op107_Nanepedanus_rufus    3.65%    failed      0.00%
 105  Op104_Dibuninae            3.65%    failed      0.33%
 106  Santinezia                 1.83%    failed      1.04%
 107  As028_Gab                  0.00%    failed      3.01%
 108  Lomanius_DNA104935         1.83%    failed      0.03%
 109  Santobius_DNA104931        1.83%    failed      0.00%
 110  As105_Phil                 0.00%    failed      3.88%
 111  As083_Cam                  0.00%    passed     15.56%
 112  Op049_Beloniscus           3.20%    passed     28.40%
 113  Bunofagea                  2.74%    failed      0.00%
 114  As031_Gab                  0.00%    failed      0.00%
 115  Fissiphallius              1.37%    failed      3.69%
 116  Larifuga                   2.74%    failed      0.17%
 117  Triaenobunus               1.37%    passed      8.13%
 118  Triaenonychidae            2.28%    failed      0.04%
 119  Equitius                   1.37%    passed     23.57%
 120  As119_Indo                 1.37%    passed     13.54%
****  TOTAL                      2.05%  79 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.013 seconds
NOTE: ModelFinder requires 19 MB RAM!
ModelFinder will test 286 DNA models (sample size: 657) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  JC            49430.810    237 99335.620    99604.861    100399.201
  2  JC+I          45835.982    238 92147.963    92420.126    93216.032
  3  JC+G4         43223.460    238 86922.920    87195.082    87990.989
  4  JC+I+G4       43020.955    239 86519.910    86795.018    87592.467
  5  JC+R2         43882.125    239 88242.250    88517.358    89314.807
  6  JC+R3         43416.702    241 87315.404    87596.474    88396.936
  7  JC+R4         43137.168    243 86760.336    87047.464    87850.843
  8  JC+R5         43002.534    245 86495.068    86788.352    87594.550
  9  JC+R6         42960.891    247 86415.782    86715.322    87524.240
 10  JC+R7         42940.836    249 86379.672    86685.568    87497.105
 11  JC+R8         42930.014    251 86362.027    86674.383    87488.436
 12  JC+R9         42923.259    253 86352.519    86671.437    87487.903
 13  JC+R10        42918.237    255 86346.474    86672.060    87490.833
 14  F81+F         48976.285    240 98432.569    98710.646    99509.613
 15  F81+F+I       45257.303    241 90996.605    91277.675    92078.137
 16  F81+F+G4      42403.193    241 85288.386    85569.456    86369.917
 17  F81+F+I+G4    42187.185    242 84858.369    85142.456    85944.389
 18  F81+F+R2      43181.345    242 86846.690    87130.777    87932.709
 19  F81+F+R3      42518.673    244 85525.345    85815.539    86620.340
 20  F81+F+R4      42222.756    246 84937.512    85233.912    86041.482
 21  F81+F+R5      42118.734    248 84733.468    85036.174    85846.414
 22  F81+F+R6      42098.250    250 84696.499    85005.612    85818.420
 23  F81+F+R7      42092.058    252 84688.115    85003.739    85819.011
 27  K2P           48701.853    238 97879.706    98151.869    98947.775
 28  K2P+I         45072.527    239 90623.054    90898.161    91695.610
 29  K2P+G4        42311.560    239 85101.119    85376.227    86173.676
 30  K2P+I+G4      42083.375    240 84646.749    84924.826    85723.793
 31  K2P+R2        43070.502    240 86621.004    86899.081    87698.048
 32  K2P+R3        42342.484    242 85168.967    85453.054    86254.987
 33  K2P+R4        42066.386    244 84620.772    84910.967    85715.767
 34  K2P+R5        42005.078    246 84502.156    84798.556    85606.126
 35  K2P+R6        41986.930    248 84469.861    84772.567    85582.806
 36  K2P+R7        41986.818    250 84473.635    84782.749    85595.556
 40  HKY+F         48154.075    241 96790.150    97071.219    97871.681
 41  HKY+F+I       44354.131    242 89192.262    89476.349    90278.281
 42  HKY+F+G4      41036.459    242 82556.919    82841.006    83642.938
 43  HKY+F+I+G4    40864.875    243 82215.750    82502.879    83306.257
 44  HKY+F+R2      42119.376    243 84724.752    85011.881    85815.260
 45  HKY+F+R3      41190.789    245 82871.579    83164.864    83971.062
 46  HKY+F+R4      40901.253    247 82296.507    82596.047    83404.965
 47  HKY+F+R5      40811.612    249 82121.224    82427.121    83238.658
 48  HKY+F+R6      40775.006    251 82052.013    82364.368    83178.422
 49  HKY+F+R7      40774.806    253 82055.613    82374.531    83190.997
 53  TNe           48543.495    239 97564.990    97840.098    98637.547
 54  TNe+I         44926.528    240 90333.056    90611.133    91410.101
 55  TNe+G4        41941.248    240 84362.496    84640.573    85439.540
 56  TNe+I+G4      41794.555    241 84071.109    84352.179    85152.641
 57  TNe+R2        42844.006    241 86170.013    86451.083    87251.545
 58  TNe+R3        42021.947    243 84529.894    84817.022    85620.401
 59  TNe+R4        41769.786    245 84029.573    84322.858    85129.056
 60  TNe+R5        41686.602    247 83867.203    84166.743    84975.661
 61  TNe+R6        41673.665    249 83845.329    84151.226    84962.763
 62  TNe+R7        41670.211    251 83842.422    84154.778    84968.831
 66  TN+F          48140.610    242 96765.219    97049.306    97851.239
 67  TN+F+I        44348.424    243 89182.848    89469.976    90273.355
 68  TN+F+G4       40979.286    243 82444.573    82731.701    83535.080
 69  TN+F+I+G4     40831.601    244 82151.202    82441.396    83246.197
 70  TN+F+R2       42089.438    244 84666.877    84957.071    85761.871
 71  TN+F+R3       41130.484    246 82752.968    83049.368    83856.939
 72  TN+F+R4       40850.318    248 82196.637    82499.343    83309.583
 73  TN+F+R5       40749.112    250 81998.224    82307.337    83120.145
 74  TN+F+R6       40726.035    252 81956.070    82271.693    83086.966
 75  TN+F+R7       40721.891    254 81951.781    82274.020    83091.653
 79  K3P           48508.013    239 97494.027    97769.135    98566.583
 80  K3P+I         44855.519    240 90191.038    90469.115    91268.083
 81  K3P+G4        42097.664    240 84675.328    84953.405    85752.372
 82  K3P+I+G4      41874.501    241 84231.001    84512.071    85312.533
 83  K3P+R2        42830.436    241 86142.872    86423.942    87224.404
 84  K3P+R3        42092.773    243 84671.547    84958.675    85762.054
 85  K3P+R4        41839.769    245 84169.538    84462.823    85269.021
 86  K3P+R5        41786.673    247 84067.347    84366.887    85175.804
 87  K3P+R6        41770.310    249 84038.620    84344.517    85156.054
 88  K3P+R7        41769.334    251 84040.667    84353.023    85167.076
 92  K3Pu+F        48102.186    242 96688.373    96972.459    97774.392
 93  K3Pu+F+I      44316.489    243 89118.978    89406.106    90209.485
 94  K3Pu+F+G4     41032.766    243 82551.532    82838.661    83642.039
 95  K3Pu+F+I+G4   40861.333    244 82210.666    82500.860    83305.660
 96  K3Pu+F+R2     42089.275    244 84666.551    84956.745    85761.546
 97  K3Pu+F+R3     41182.829    246 82857.657    83154.057    83961.627
 98  K3Pu+F+R4     40894.786    248 82285.572    82588.278    83398.518
 99  K3Pu+F+R5     40801.448    250 82102.896    82412.009    83224.817
100  K3Pu+F+R6     40772.201    252 82048.403    82364.026    83179.299
101  K3Pu+F+R7     40770.003    254 82048.007    82370.245    83187.878
105  TPM2+F        47804.816    242 96093.632    96377.719    97179.651
106  TPM2+F+I      44133.319    243 88752.638    89039.766    89843.145
107  TPM2+F+G4     41027.727    243 82541.453    82828.581    83631.960
108  TPM2+F+I+G4   40853.210    244 82194.421    82484.615    83289.416
109  TPM2+F+R2     41994.952    244 84477.904    84768.098    85572.899
110  TPM2+F+R3     41162.201    246 82816.402    83112.802    83920.372
111  TPM2+F+R4     40875.671    248 82247.341    82550.047    83360.287
112  TPM2+F+R5     40791.562    250 82083.125    82392.238    83205.046
113  TPM2+F+R6     40766.279    252 82036.557    82352.181    83167.454
114  TPM2+F+R7     40763.283    254 82034.566    82356.805    83174.438
118  TPM2u+F       47804.815    242 96093.630    96377.717    97179.649
119  TPM2u+F+I     44133.319    243 88752.638    89039.766    89843.145
120  TPM2u+F+G4    41027.726    243 82541.452    82828.581    83631.960
121  TPM2u+F+I+G4  40853.221    244 82194.441    82484.635    83289.436
122  TPM2u+F+R2    41994.887    244 84477.774    84767.968    85572.768
123  TPM2u+F+R3    41162.200    246 82816.400    83112.800    83920.370
124  TPM2u+F+R4    40875.677    248 82247.354    82550.060    83360.299
125  TPM2u+F+R5    40790.506    250 82081.011    82390.124    83202.932
126  TPM2u+F+R6    40766.368    252 82036.737    82352.361    83167.633
127  TPM2u+F+R7    40763.159    254 82034.318    82356.557    83174.190
131  TPM3+F        48116.966    242 96717.932    97002.019    97803.951
132  TPM3+F+I      44317.132    243 89120.264    89407.392    90210.771
133  TPM3+F+G4     41025.251    243 82536.503    82823.631    83627.010
134  TPM3+F+I+G4   40853.019    244 82194.037    82484.231    83289.032
135  TPM3+F+R2     42090.926    244 84669.851    84960.046    85764.846
136  TPM3+F+R3     41175.451    246 82842.901    83139.301    83946.871
137  TPM3+F+R4     40887.552    248 82271.103    82573.809    83384.049
138  TPM3+F+R5     40790.574    250 82081.148    82390.262    83203.069
139  TPM3+F+R6     40764.372    252 82032.745    82348.368    83163.641
140  TPM3+F+R7     40762.051    254 82032.103    82354.342    83171.975
144  TPM3u+F       48116.966    242 96717.932    97002.019    97803.952
145  TPM3u+F+I     44317.132    243 89120.264    89407.392    90210.771
146  TPM3u+F+G4    41025.252    243 82536.504    82823.632    83627.011
147  TPM3u+F+I+G4  40853.022    244 82194.044    82484.238    83289.039
148  TPM3u+F+R2    42090.905    244 84669.810    84960.004    85764.805
149  TPM3u+F+R3    41175.437    246 82842.873    83139.273    83946.844
150  TPM3u+F+R4    40887.559    248 82271.119    82573.825    83384.064
151  TPM3u+F+R5    40790.619    250 82081.237    82390.350    83203.158
152  TPM3u+F+R6    40764.358    252 82032.715    82348.339    83163.612
153  TPM3u+F+R7    40761.942    254 82031.884    82354.123    83171.756
157  TIMe          48348.229    240 97176.457    97454.534    98253.501
158  TIMe+I        44708.287    241 89898.574    90179.644    90980.106
159  TIMe+G4       41753.731    241 83989.463    84270.533    85070.995
160  TIMe+I+G4     41589.336    242 83662.672    83946.758    84748.691
161  TIMe+R2       42607.854    242 85699.708    85983.795    86785.728
162  TIMe+R3       41798.420    244 84084.840    84375.034    85179.834
163  TIMe+R4       41552.957    246 83597.915    83894.315    84701.885
164  TIMe+R5       41488.941    248 83473.882    83776.588    84586.828
165  TIMe+R6       41478.226    250 83456.452    83765.566    84578.373
166  TIMe+R7       41469.371    252 83442.742    83758.365    84573.638
167  TIMe+R8       41467.731    254 83443.462    83765.701    84583.334
170  TIM+F         48088.698    243 96663.396    96950.524    97753.903
171  TIM+F+I       44310.698    244 89109.396    89399.590    90204.391
172  TIM+F+G4      40975.629    244 82439.258    82729.452    83534.253
173  TIM+F+I+G4    40826.917    245 82143.833    82437.118    83243.316
174  TIM+F+R2      42063.201    245 84616.403    84909.688    85715.885
175  TIM+F+R3      41124.691    247 82743.381    83042.922    83851.839
176  TIM+F+R4      40843.976    249 82185.951    82491.848    83303.384
177  TIM+F+R5      40742.372    251 81986.744    82299.100    83113.153
178  TIM+F+R6      40723.647    253 81953.293    82272.212    83088.677
179  TIM+F+R7      40711.106    255 81932.212    82257.798    83076.572
180  TIM+F+R8      40709.341    257 81932.681    82265.042    83086.016
183  TIM2e         47774.728    240 96029.455    96307.532    97106.500
184  TIM2e+I       44307.426    241 89096.851    89377.921    90178.383
185  TIM2e+G4      41628.312    241 83738.625    84019.695    84820.157
186  TIM2e+I+G4    41468.919    242 83421.838    83705.925    84507.857
187  TIM2e+R2      42370.407    242 85224.814    85508.901    86310.833
188  TIM2e+R3      41648.238    244 83784.477    84074.671    84879.472
189  TIM2e+R4      41424.875    246 83341.751    83638.151    84445.721
190  TIM2e+R5      41358.066    248 83212.133    83514.839    84325.078
191  TIM2e+R6      41350.018    250 83200.037    83509.150    84321.958
192  TIM2e+R7      41343.401    252 83190.801    83506.425    84321.698
193  TIM2e+R8      41343.277    254 83194.555    83516.794    84334.427
196  TIM2+F        47795.601    243 96077.201    96364.329    97167.708
197  TIM2+F+I      44128.982    244 88745.964    89036.158    89840.959
198  TIM2+F+G4     40966.687    244 82421.374    82711.568    83516.369
199  TIM2+F+I+G4   40808.285    245 82106.570    82399.855    83206.053
200  TIM2+F+R2     41970.605    245 84431.210    84724.495    85530.693
201  TIM2+F+R3     41100.854    247 82695.708    82995.248    83804.165
202  TIM2+F+R4     40820.061    249 82138.123    82444.019    83255.556
203  TIM2+F+R5     40725.780    251 81953.561    82265.916    83079.970
204  TIM2+F+R6     40706.246    253 81918.492    82237.410    83053.876
205  TIM2+F+R7     40694.307    255 81898.613    82224.199    83042.972
206  TIM2+F+R8     40693.980    257 81901.959    82234.320    83055.294
209  TIM3e         48366.863    240 97213.726    97491.803    98290.770
210  TIM3e+I       44767.340    241 90016.680    90297.750    91098.212
211  TIM3e+G4      41848.915    241 84179.831    84460.901    85261.363
212  TIM3e+I+G4    41692.869    242 83869.738    84153.825    84955.758
213  TIM3e+R2      42693.325    242 85870.650    86154.737    86956.670
214  TIM3e+R3      41913.941    244 84315.881    84606.075    85410.876
215  TIM3e+R4      41659.505    246 83811.011    84107.411    84914.981
216  TIM3e+R5      41587.646    248 83671.292    83973.998    84784.237
217  TIM3e+R6      41573.288    250 83646.576    83955.689    84768.497
218  TIM3e+R7      41568.873    252 83641.747    83957.370    84772.643
222  TIM3+F        48102.935    243 96691.870    96978.998    97782.377
223  TIM3+F+I      44311.100    244 89110.201    89400.395    90205.196
224  TIM3+F+G4     40970.843    244 82429.687    82719.881    83524.681
225  TIM3+F+I+G4   40822.768    245 82135.535    82428.820    83235.018
226  TIM3+F+R2     42066.752    245 84623.505    84916.789    85722.987
227  TIM3+F+R3     41121.684    247 82737.368    83036.908    83845.826
228  TIM3+F+R4     40840.754    249 82179.508    82485.405    83296.942
229  TIM3+F+R5     40738.114    251 81978.228    82290.583    83104.636
230  TIM3+F+R6     40719.989    253 81945.977    82264.895    83081.361
231  TIM3+F+R7     40707.577    255 81925.155    82250.741    83069.514
232  TIM3+F+R8     40704.198    257 81922.396    82254.757    83075.730
235  TVMe          47728.010    241 95938.021    96219.091    97019.553
236  TVMe+I        44257.550    242 88999.100    89283.187    90085.119
237  TVMe+G4       41780.060    242 84044.119    84328.206    85130.139
238  TVMe+I+G4     41576.608    243 83639.215    83926.343    84729.722
239  TVMe+R2       42399.850    243 85285.700    85572.828    86376.207
240  TVMe+R3       41740.465    245 83970.930    84264.215    85070.413
241  TVMe+R4       41520.419    247 83534.838    83834.378    84643.296
242  TVMe+R5       41474.085    249 83446.169    83752.066    84563.603
243  TVMe+R6       41467.396    251 83436.792    83749.148    84563.201
244  TVMe+R7       41464.363    253 83434.725    83753.643    84570.109
248  TVM+F         47767.625    244 96023.249    96313.443    97118.244
249  TVM+F+I       44095.880    245 88681.759    88975.044    89781.242
250  TVM+F+G4      41016.563    245 82523.126    82816.411    83622.609
251  TVM+F+I+G4    40841.869    246 82175.739    82472.139    83279.709
252  TVM+F+R2      41977.985    246 84447.969    84744.369    85551.939
253  TVM+F+R3      41149.681    248 82795.361    83098.067    83908.307
254  TVM+F+R4      40865.103    250 82230.206    82539.319    83352.127
255  TVM+F+R5      40778.259    252 82060.518    82376.142    83191.414
256  TVM+F+R6      40760.478    254 82028.956    82351.195    83168.828
257  TVM+F+R7      40750.074    256 82012.148    82341.108    83160.995
258  TVM+F+R8      40744.321    258 82004.641    82340.430    83162.464
261  SYM           47586.146    242 95656.291    95940.378    96742.311
262  SYM+I         44120.376    243 88726.752    89013.880    89817.259
263  SYM+G4        41456.062    243 83398.124    83685.252    84488.631
264  SYM+I+G4      41291.931    244 83071.862    83362.056    84166.857
265  SYM+R2        42180.010    244 84848.020    85138.214    85943.015
266  SYM+R3        41470.142    246 83432.284    83728.684    84536.255
267  SYM+R4        41242.914    248 82981.828    83284.534    84094.774
268  SYM+R5        41188.225    250 82876.451    83185.564    83998.372
269  SYM+R6        41178.279    252 82860.557    83176.181    83991.454
270  SYM+R7        41176.207    254 82860.415    83182.654    84000.287
274  GTR+F         47757.983    245 96005.965    96299.250    97105.448
275  GTR+F+I       44091.263    246 88674.526    88970.926    89778.496
276  GTR+F+G4      40958.489    246 82408.977    82705.377    83512.947
277  GTR+F+I+G4    40800.305    247 82094.609    82394.150    83203.067
278  GTR+F+R2      41954.747    247 84403.495    84703.035    85511.952
279  GTR+F+R3      41092.792    249 82683.584    82989.481    83801.017
280  GTR+F+R4      40812.176    251 82126.351    82438.707    83252.760
281  GTR+F+R5      40717.005    253 81940.010    82258.928    83075.394
282  GTR+F+R6      40701.283    255 81912.567    82238.153    83056.926
283  GTR+F+R7      40693.875    257 81901.751    82234.112    83055.086
284  GTR+F+R8      40687.716    259 81893.433    82232.677    83055.743
Akaike Information Criterion:           GTR+F+R8
Corrected Akaike Information Criterion: TIM2+F+R7
Bayesian Information Criterion:         TIM2+F+R7
Best-fit model: TIM2+F+R7 chosen according to BIC

All model information printed to COI-aligned-mafft.fasta.model.gz
CPU time for ModelFinder: 291.604 seconds (0h:4m:51s)
Wall-clock time for ModelFinder: 292.054 seconds (0h:4m:52s)

NOTE: 13 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -42600.195
2. Current log-likelihood: -40995.120
3. Current log-likelihood: -40897.318
4. Current log-likelihood: -40869.950
5. Current log-likelihood: -40863.028
6. Current log-likelihood: -40858.967
7. Current log-likelihood: -40858.076
8. Current log-likelihood: -40856.567
9. Current log-likelihood: -40856.304
10. Current log-likelihood: -40855.187
11. Current log-likelihood: -40855.082
12. Current log-likelihood: -40854.532
Optimal log-likelihood: -40854.485
Rate parameters:  A-C: 1.51558  A-G: 3.62533  A-T: 1.51558  C-G: 1.00000  C-T: 5.65204  G-T: 1.00000
Base frequencies:  A: 0.270  C: 0.215  G: 0.147  T: 0.369
Site proportion and rates:  (0.395,0.050) (0.033,0.307) (0.047,0.309) (0.076,0.653) (0.177,1.097) (0.125,1.824) (0.146,3.307)
Parameters optimization took 12 rounds (4.592 sec)
Computing ML distances based on estimated model parameters... 0.149 sec
Computing BIONJ tree...
0.005 seconds
Log-likelihood of BIONJ tree: -40669.824
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 1.006 second
Computing log-likelihood of 98 initial trees ... 3.712 seconds
Current best score: -40669.824

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -40522.835
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 7: -40491.891
Iteration 10 / LogL: -40536.655 / Time: 0h:0m:21s
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 11: -40468.593
Iteration 20 / LogL: -40521.311 / Time: 0h:0m:28s
Finish initializing candidate tree set (20)
Current best tree score: -40468.593 / CPU time: 23.230
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 28: -40441.528
Iteration 30 / LogL: -40474.261 / Time: 0h:0m:35s (0h:2m:0s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 35: -40430.654
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 36: -40425.364
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 39: -40412.855
Iteration 40 / LogL: -40482.965 / Time: 0h:0m:45s (0h:1m:54s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 41: -40407.081
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 43: -40387.846
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 47: -40380.270
Iteration 50 / LogL: -40409.190 / Time: 0h:0m:54s (0h:1m:48s left)
Iteration 60 / LogL: -40391.935 / Time: 0h:0m:58s (0h:1m:26s left)
Iteration 70 / LogL: -40382.404 / Time: 0h:1m:3s (0h:1m:10s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 75: -40377.953
Iteration 80 / LogL: -40390.752 / Time: 0h:1m:7s (0h:1m:21s left)
Iteration 90 / LogL: -40381.491 / Time: 0h:1m:11s (0h:1m:8s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 95: -40375.214
Iteration 100 / LogL: -40376.889 / Time: 0h:1m:16s (0h:1m:13s left)
Iteration 110 / LogL: -40396.062 / Time: 0h:1m:20s (0h:1m:2s left)
Iteration 120 / LogL: -40386.573 / Time: 0h:1m:24s (0h:0m:53s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 128: -40374.740
Iteration 130 / LogL: -40406.888 / Time: 0h:1m:29s (0h:1m:7s left)
Iteration 140 / LogL: -40385.533 / Time: 0h:1m:33s (0h:0m:59s left)
Iteration 150 / LogL: -40388.964 / Time: 0h:1m:37s (0h:0m:50s left)
Iteration 160 / LogL: -40387.623 / Time: 0h:1m:41s (0h:0m:43s left)
Iteration 170 / LogL: -40388.092 / Time: 0h:1m:45s (0h:0m:36s left)
Iteration 180 / LogL: -40382.199 / Time: 0h:1m:50s (0h:0m:29s left)
Iteration 190 / LogL: -40382.393 / Time: 0h:1m:54s (0h:0m:22s left)
Iteration 200 / LogL: -40377.961 / Time: 0h:1m:58s (0h:0m:16s left)
Iteration 210 / LogL: -40392.158 / Time: 0h:2m:3s (0h:0m:10s left)
Iteration 220 / LogL: -40386.158 / Time: 0h:2m:7s (0h:0m:4s left)
TREE SEARCH COMPLETED AFTER 229 ITERATIONS / Time: 0h:2m:11s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -40374.740
Optimal log-likelihood: -40374.732
Rate parameters:  A-C: 1.47496  A-G: 3.68704  A-T: 1.47496  C-G: 1.00000  C-T: 6.83298  G-T: 1.00000
Base frequencies:  A: 0.270  C: 0.215  G: 0.147  T: 0.369
Site proportion and rates:  (0.369,0.019) (0.072,0.176) (0.052,0.331) (0.082,0.645) (0.187,1.077) (0.122,2.006) (0.116,3.999)
Parameters optimization took 1 rounds (0.333 sec)
BEST SCORE FOUND : -40374.732
Total tree length: 35.871

Total number of iterations: 229
CPU time used for tree search: 126.277 sec (0h:2m:6s)
Wall-clock time used for tree search: 126.450 sec (0h:2m:6s)
Total CPU time used: 131.488 sec (0h:2m:11s)
Total wall-clock time used: 131.677 sec (0h:2m:11s)

Analysis results written to: 
  IQ-TREE report:                COI-aligned-mafft.fasta.iqtree
  Maximum-likelihood tree:       COI-aligned-mafft.fasta.treefile
  Likelihood distances:          COI-aligned-mafft.fasta.mldist
  Screen log file:               COI-aligned-mafft.fasta.log
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputss/COI-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies/Gene_Trees

#### Histone H3 Gene Tree
<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s H3-aligned-mafft.fasta>
-running H3 Mafft Alignment
Output:
Reading alignment file H3-aligned-mafft.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
Alignment has 151 sequences with 327 columns, 219 distinct patterns
137 parsimony-informative, 20 singleton sites, 169 constant sites
                         Gap/Ambiguity  Composition  p-value
   1  Agoristenidae_DNA105839    0.00%    passed      7.91%
   2  Assamiidae_DNA104857       0.00%    passed     66.46%
   3  As109_Thai                 0.00%    passed     69.83%
   4  As120_Thai                 0.00%    passed     66.22%
   5  As125_Thai                 0.31%    passed     92.53%
   6  Assamiidae_DNA104859       0.00%    passed     51.73%
   7  As110_Laos                 2.45%    passed     47.71%
   8  As121_Laos                 0.61%    passed     51.54%
   9  As132_Mala                 1.53%    passed     80.30%
  10  As108_Laos                 0.61%    passed     76.47%
  11  As114_Thai                 0.61%    passed     99.33%
  12  Chilon                     0.00%    passed     63.54%
  13  As034_Cam                  9.48%    passed     69.93%
  14  As027_Gab                 49.24%    passed     74.41%
  15  As085_Cam                  9.48%    passed     62.39%
  16  As084_Cam                 12.84%    passed     63.07%
  17  As056_Gab                  8.56%    passed     95.35%
  18  As082_Cam                 10.09%    passed     92.92%
  19  As032_Gab                 64.83%    passed     89.49%
  20  As039_Gab                 12.23%    passed     72.26%
  21  As021_Gab                  9.48%    passed     98.63%
  22  As072_Lib                 11.01%    passed     65.04%
  23  As028_Gab                 19.27%    passed     79.17%
  24  As045_Cam                  8.56%    passed     81.96%
  25  As070_Lib                 10.40%    passed     58.74%
  26  As083_Cam                  9.79%    passed     92.12%
  27  As014_Gab                  8.26%    passed     81.31%
  28  As059_Gab                  9.48%    passed     92.79%
  29  Paraselenca                0.00%    passed     99.66%
  30  As040_Gab                  7.65%    passed     97.18%
  31  As057_Gab                 19.27%    passed     90.33%
  32  As018_Gab                 25.69%    passed     55.61%
  33  As050_Cam                  7.03%    passed     70.22%
  34  As017_Gab                 10.40%    passed     75.71%
  35  As116_Thai                 0.31%    passed     94.83%
  36  As117_Thai                 2.14%    passed     98.58%
  37  Neopygoplus                0.00%    passed     75.72%
  38  As127_Laos                 0.92%    passed     93.16%
  39  As133_Indo                 0.92%    passed     99.46%
  40  As124_Indo                 9.79%    passed     89.98%
  41  As080_Indo                 9.79%    passed     96.84%
  42  As081_Indo                 8.56%    passed     83.01%
  43  As130_Viet                 0.00%    passed     97.50%
  44  As134_Viet                 0.31%    passed     97.15%
  45  As085_WAus                 9.79%    passed     54.67%
  46  As087_WAus                 7.95%    passed     77.11%
  47  As101_NAus                 9.17%    passed     89.89%
  48  As102_NAus                 8.87%    passed     80.03%
  49  As103_NAus                10.40%    passed     58.37%
  50  As089_EAus                 8.56%    passed     88.87%
  51  As092_EAus                55.35%    passed     11.30%
  52  As104_EAus                 8.56%    passed     82.95%
  53  As094_PNG                 22.63%    passed     37.77%
  54  As095_PNG                 26.61%    passed     25.67%
  55  As096_PNG                 12.84%    passed     55.64%
  56  As097_PNG                  9.79%    passed     93.66%
  57  Epedanidae_DNA104062       0.00%    passed     97.21%
  58  As119_Indo                 0.31%    passed     68.69%
  59  Epedanidae_DNA104066       0.00%    passed     85.68%
  60  Tithaeus                   0.00%    passed     99.77%
  61  As131_Mala                 0.31%    passed     96.63%
  62  As115_Thai                 0.31%    passed     99.39%
  63  As118_Thai                 0.31%    passed     94.39%
  64  As122_Laos                 2.45%    passed     98.90%
  65  As129_Laos                 0.92%    passed     99.05%
  66  As136_Viet                 0.61%    passed     98.67%
  67  As140_Viet                 0.61%    passed     89.14%
  68  As126_Laos                 0.61%    passed     99.71%
  69  As123_Viet                 0.31%    passed     91.69%
  70  As106_Thai                 0.31%    passed     97.97%
  71  Conomma                    0.00%    passed     51.70%
  72  As020_Gab                 59.02%    passed     83.85%
  73  As030_Gab                 26.30%    passed     66.54%
  74  As041_Gab                 12.54%    passed     64.00%
  75  As038_Gab                  9.48%    passed     99.01%
  76  As026_Gab                  8.26%    passed     83.67%
  77  Stygnopsidae_DNA103882     0.00%    passed     93.56%
  78  Stygnopsidae_DNA104855     0.00%    passed     88.23%
  79  Stygnopsidae_DNA104856     0.00%    passed     89.01%
  80  Synthetonychia             0.00%    passed     22.09%
  81  Fumontana                  0.00%    passed     44.63%
  82  Fijicolana                 0.00%    passed     86.31%
  83  Pellobunus                 0.00%    passed     94.47%
  84  As061_Lib                  8.26%    passed     88.06%
  85  As069_Lib                 18.65%    passed     72.22%
  86  Cynortula                  0.00%    passed     65.38%
  87  Metalibitia                0.31%    passed     81.19%
  88  Zygopachylus               1.53%    passed     89.19%
  89  Glysterus                  3.06%    passed     97.49%
  90  Heterocranaus              0.00%    passed     96.28%
  91  Santinezia                 0.31%    passed     64.19%
  92  Stygnoplus                 0.00%    passed     54.90%
  93  Epedanidae_DNA104861       0.00%    passed     94.88%
  94  Epedanidae_DNA104862       0.00%    passed     80.50%
  95  Pseudoepedanus             0.00%    passed     82.81%
  96  Stygnomma                  0.00%    passed     89.80%
  97  Guasinia                   0.00%    passed     48.47%
  98  Kimula                    27.83%    passed     78.38%
  99  Minuella                   2.45%    passed     72.07%
 100  IC_DNA102668               0.31%    passed     37.98%
 101  IC_DNA102669               2.14%    passed     29.43%
 102  Equitius                   0.00%    passed     98.79%
 103  Triaenobunus               0.00%    passed     77.81%
 104  Triaenonychidae            0.00%    passed     99.91%
 105  Peltonychia                0.00%    passed     87.77%
 106  Theromaster                0.00%    passed     37.89%
 107  As098_Phil                 9.79%    passed     50.25%
 108  As105_Phil                14.07%    passed     89.38%
 109  Hesperonemastoma           0.00%    passed     40.37%
 110  Pantopsalis                2.75%    passed     50.46%
 111  Zalmoxida                  1.53%    passed     76.69%
 112  Trojanella                 0.00%    passed     52.10%
 113  As135_Viet                 0.31%    passed     97.15%
 114  Trionyxella                0.92%    passed     37.45%
 115  As111_Phil                 0.92%    passed     92.34%
 116  Biantidae_DNA105668        0.31%    passed     97.97%
 117  Baculigerus                0.00%    passed     37.57%
 118  Escadabius                 0.00%    passed     83.23%
 119  Ethobunus                  0.00%    passed     57.11%
 120  Zalmoxis                   0.00%    passed     58.02%
 121  Icaleptes                  0.00%    passed     89.80%
 122  Icaleptes_DNA104053        0.00%    passed     89.99%
 123  Icaleptes_DNA104056        0.00%    passed     86.62%
 124  Fissiphallius              0.00%    passed     77.92%
 125  Fissiphallius_DNA104055    0.00%    passed     78.06%
 126  Fissiphallius_DNA104057    0.61%    passed     76.63%
 127  Larifuga                   0.00%    passed     44.75%
 128  Metabiantes_DNA100704      0.61%    passed     65.69%
 129  Samoidae                   0.00%    passed     72.18%
 130  Haasus_sp                  0.00%    passed     87.96%
 131  Haasus_judaeus             0.00%    passed     60.07%
 132  Holoscotolemon             0.00%    passed     96.47%
 133  Caenoncopus                0.00%    passed     55.28%
 134  Gnomulus                   0.00%    passed     31.42%
 135  Palaeoncopus               0.31%    passed     86.50%
 136  Sandokan                   0.61%    passed     96.22%
 137  Martensiellus              0.00%    passed     53.71%
 138  Arulla_DNA102666           0.61%    passed     64.99%
 139  Hoplobunus                 0.00%    passed     82.87%
 140  Karos                      0.31%    passed     74.24%
 141  Troglosiro                 0.61%    failed      1.19%
 142  Rostromontia               0.00%    passed     60.54%
 143  Dongmoa                    1.22%    passed     21.63%
 144  Lomanius_DNA104934         0.00%    passed     57.36%
 145  Lomanius_DNA104935         1.22%    passed     69.54%
 146  IC_DNA103572               0.61%    passed     95.70%
 147  Megapachylus               2.14%    passed      9.47%
 148  Protolophus                0.00%    failed      0.15%
 149  Bunofagea                  0.61%    passed     23.00%
 150  Scotolemon_lespesi         0.00%    passed     33.14%
 151  Goniosoma                  0.00%    failed      1.85%
WARNING: 3 sequences contain more than 50% gaps/ambiguity
****  TOTAL                      5.23%  3 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.010 seconds
NOTE: ModelFinder requires 10 MB RAM!
ModelFinder will test 286 DNA models (sample size: 327) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  JC            13997.923    299 28593.845    35238.290    29727.043
  2  JC+I          12174.922    300 24949.845    31895.999    26086.833
  3  JC+G4         11544.660    300 23689.320    30635.474    24826.308
  4  JC+I+G4       11428.035    301 23458.071    30730.231    24598.849
  5  JC+R2         11893.567    301 24389.133    31661.293    25529.911
  6  JC+R3         11585.840    303 23777.681    31787.420    24926.039
  7  JC+R4         11496.788    305 23603.576    32492.147    24759.513
  8  JC+R5         11474.414    307 23562.827    33516.090    24726.345
  9  JC+R6         11457.336    309 23532.672    34802.084    24703.770
 10  JC+R7         11446.034    311 23514.068    36451.668    24692.746
 11  JC+R8         11437.911    313 23501.823    38622.131    24688.080
 12  JC+R9         11433.776    315 23497.552    41595.734    24691.390
 14  F81+F         14151.555    302 28907.110    36532.610    30051.678
 15  F81+F+I       12267.535    303 25141.071    33150.810    26289.428
 16  F81+F+G4      11638.545    303 23883.090    31892.829    25031.448
 17  F81+F+I+G4    11521.352    304 23650.704    32079.795    24802.852
 18  F81+F+R2      11949.311    304 24506.621    32935.712    25658.769
 19  F81+F+R3      11632.108    306 23876.217    33270.417    25035.945
 20  F81+F+R4      11549.497    308 23714.994    34289.661    24882.302
 21  F81+F+R5      11533.682    310 23687.364    35738.614    24862.252
 22  F81+F+R6      11524.411    312 23672.822    37623.679    24855.289
 23  F81+F+R7      11519.311    314 23666.623    40151.623    24856.670
 27  K2P           13766.203    300 28132.407    35078.560    29269.395
 28  K2P+I         11932.146    301 24466.291    31738.451    25607.069
 29  K2P+G4        11281.491    301 23164.982    30437.142    24305.760
 30  K2P+I+G4      11158.110    302 22920.220    30545.720    24064.788
 31  K2P+R2        11555.923    302 23715.846    31341.346    24860.414
 32  K2P+R3        11239.486    304 23086.973    31516.063    24239.120
 33  K2P+R4        11165.825    306 22943.649    32337.849    24103.377
 34  K2P+R5        11155.226    308 22926.451    33501.118    24093.759
 35  K2P+R6        11152.347    310 22924.693    34975.943    24099.581
 40  HKY+F         13880.861    303 28367.721    36377.460    29516.079
 41  HKY+F+I       11933.571    304 24475.141    32904.232    25627.289
 42  HKY+F+G4      11279.075    304 23166.149    31595.240    24318.297
 43  HKY+F+I+G4    11149.196    305 22908.391    31796.963    24064.329
 44  HKY+F+R2      11552.680    305 23715.360    32603.931    24871.298
 45  HKY+F+R3      11215.963    307 23045.925    32999.188    24209.443
 46  HKY+F+R4      11150.868    309 22919.736    34189.148    24090.834
 47  HKY+F+R5      11140.707    311 22903.414    35841.014    24082.092
 48  HKY+F+R6      11136.452    313 22898.904    38019.211    24085.161
 53  TNe           13708.660    301 28019.320    35291.480    29160.098
 54  TNe+I         11929.318    302 24462.636    32088.136    25607.204
 55  TNe+G4        11280.964    302 23165.928    30791.428    24310.496
 56  TNe+I+G4      11157.723    303 22921.446    30931.185    24069.804
 57  TNe+R2        11508.884    303 23623.769    31633.508    24772.127
 58  TNe+R3        11221.512    305 23053.023    31941.595    24208.961
 59  TNe+R4        11155.695    307 22925.390    32878.653    24088.907
 60  TNe+R5        11150.839    309 22919.679    34189.090    24090.776
 66  TN+F          13790.609    304 28189.218    36618.309    29341.366
 67  TN+F+I        11910.638    305 24431.277    33319.848    25587.214
 68  TN+F+G4       11270.497    305 23150.994    32039.566    24306.932
 69  TN+F+I+G4     11139.282    306 22890.565    32284.765    24050.293
 70  TN+F+R2       11493.389    306 23598.778    32992.978    24758.506
 71  TN+F+R3       11199.821    308 23015.642    33590.308    24182.949
 72  TN+F+R4       11137.793    310 22895.586    34946.836    24070.474
 73  TN+F+R5       11130.037    312 22884.074    36834.931    24066.541
 74  TN+F+R6       11126.123    314 22880.247    39365.247    24070.294
 79  K3P           13765.801    301 28133.602    35405.762    29274.380
 80  K3P+I         11931.616    302 24467.233    32092.733    25611.801
 81  K3P+G4        11281.362    302 23166.723    30792.223    24311.291
 82  K3P+I+G4      11157.990    303 22921.980    30931.719    24070.338
 83  K3P+R2        11494.333    303 23594.666    31604.405    24743.024
 84  K3P+R3        11220.260    305 23050.519    31939.091    24206.457
 85  K3P+R4        11155.789    307 22925.579    32878.842    24089.096
 86  K3P+R5        11150.018    309 22918.035    34187.447    24089.133
 92  K3Pu+F        13880.794    304 28369.588    36798.679    29521.736
 93  K3Pu+F+I      11932.379    305 24474.757    33363.329    25630.695
 94  K3Pu+F+G4     11276.277    305 23162.555    32051.126    24318.493
 95  K3Pu+F+I+G4   11145.543    306 22903.086    32297.286    24062.814
 96  K3Pu+F+R2     11494.941    306 23601.882    32996.082    24761.610
 97  K3Pu+F+R3     11203.464    308 23022.929    33597.595    24190.236
 98  K3Pu+F+R4     11146.153    310 22912.307    34963.557    24087.195
 99  K3Pu+F+R5     11134.333    312 22892.665    36843.523    24075.133
100  K3Pu+F+R6     11130.651    314 22889.301    39374.301    24079.349
105  TPM2+F        13874.818    304 28357.636    36786.727    29509.784
106  TPM2+F+I      11916.042    305 24442.083    33330.655    25598.021
107  TPM2+F+G4     11250.439    305 23110.878    31999.450    24266.816
108  TPM2+F+I+G4   11115.970    306 22843.940    32238.140    24003.668
109  TPM2+F+R2     11463.570    306 23539.140    32933.340    24698.868
110  TPM2+F+R3     11176.419    308 22968.838    33543.504    24136.145
111  TPM2+F+R4     11117.994    310 22855.987    34907.237    24030.875
112  TPM2+F+R5     11105.420    312 22834.840    36785.697    24017.307
113  TPM2+F+R6     11102.394    314 22832.788    39317.788    24022.836
118  TPM2u+F       13874.816    304 28357.632    36786.723    29509.780
119  TPM2u+F+I     11916.040    305 24442.080    33330.651    25598.018
120  TPM2u+F+G4    11250.436    305 23110.872    31999.444    24266.810
121  TPM2u+F+I+G4  11115.970    306 22843.940    32238.140    24003.668
122  TPM2u+F+R2    11456.980    306 23525.959    32920.159    24685.687
123  TPM2u+F+R3    11176.199    308 22968.399    33543.065    24135.706
124  TPM2u+F+R4    11117.715    310 22855.430    34906.680    24030.318
125  TPM2u+F+R5    11105.167    312 22834.334    36785.191    24016.801
126  TPM2u+F+R6    11102.173    314 22832.346    39317.346    24022.393
131  TPM3+F        13870.883    304 28349.765    36778.856    29501.913
132  TPM3+F+I      11932.529    305 24475.059    33363.630    25630.997
133  TPM3+F+G4     11279.001    305 23168.002    32056.573    24323.939
134  TPM3+F+I+G4   11148.499    306 22908.997    32303.197    24068.725
135  TPM3+F+R2     11479.305    306 23570.611    32964.811    24730.339
136  TPM3+F+R3     11205.558    308 23027.117    33601.784    24194.425
137  TPM3+F+R4     11148.284    310 22916.567    34967.817    24091.455
138  TPM3+F+R5     11135.606    312 22895.212    36846.069    24077.680
139  TPM3+F+R6     11132.041    314 22892.082    39377.082    24082.129
144  TPM3u+F       13870.879    304 28349.757    36778.848    29501.905
145  TPM3u+F+I     11932.526    305 24475.053    33363.624    25630.990
146  TPM3u+F+G4    11279.010    305 23168.020    32056.592    24323.958
147  TPM3u+F+I+G4  11148.494    306 22908.987    32303.187    24068.715
148  TPM3u+F+R2    11475.779    306 23563.558    32957.758    24723.286
149  TPM3u+F+R3    11205.465    308 23026.930    33601.596    24194.237
150  TPM3u+F+R4    11148.096    310 22916.192    34967.442    24091.079
151  TPM3u+F+R5    11135.532    312 22895.063    36845.920    24077.531
152  TPM3u+F+R6    11131.932    314 22891.865    39376.865    24081.912
157  TIMe          13708.295    302 28020.590    35646.090    29165.158
158  TIMe+I        11928.863    303 24463.726    32473.465    25612.084
159  TIMe+G4       11280.844    303 23167.688    31177.427    24316.046
160  TIMe+I+G4     11157.602    304 22923.205    31352.295    24075.352
161  TIMe+R2       11473.781    304 23555.562    31984.653    24707.710
162  TIMe+R3       11219.524    306 23051.048    32445.248    24210.775
163  TIMe+R4       11154.207    308 22924.413    33499.080    24091.721
164  TIMe+R5       11149.114    310 22918.228    34969.478    24093.116
170  TIM+F         13790.546    305 28191.091    37079.662    29347.029
171  TIM+F+I       11909.442    306 24430.884    33825.084    25590.612
172  TIM+F+G4      11267.822    306 23147.644    32541.844    24307.372
173  TIM+F+I+G4    11135.536    307 22885.072    32838.335    24048.590
174  TIM+F+R2      11452.817    307 23519.634    33472.897    24683.151
175  TIM+F+R3      11194.368    309 23006.735    34276.147    24177.833
176  TIM+F+R4      11132.529    311 22887.058    35824.658    24065.736
177  TIM+F+R5      11124.426    313 22874.852    37995.160    24061.110
178  TIM+F+R6      11120.768    315 22871.536    40969.718    24065.374
183  TIM2e         13706.734    302 28017.469    35642.969    29162.037
184  TIM2e+I       11922.065    303 24450.130    32459.869    25598.488
185  TIM2e+G4      11267.856    303 23141.711    31151.450    24290.069
186  TIM2e+I+G4    11143.917    304 22895.834    31324.925    24047.982
187  TIM2e+R2      11454.683    304 23517.366    31946.456    24669.513
188  TIM2e+R3      11205.122    306 23022.243    32416.443    24181.971
189  TIM2e+R4      11140.847    308 22897.694    33472.361    24065.002
190  TIM2e+R5      11135.642    310 22891.284    34942.534    24066.172
196  TIM2+F        13784.509    305 28179.017    37067.589    29334.955
197  TIM2+F+I      11893.623    306 24399.246    33793.446    25558.974
198  TIM2+F+G4     11242.469    306 23096.939    32491.139    24256.666
199  TIM2+F+I+G4   11108.442    307 22830.883    32784.147    23994.401
200  TIM2+F+R2     11418.845    307 23451.690    33404.953    24615.208
201  TIM2+F+R3     11168.005    309 22954.010    34223.422    24125.108
202  TIM2+F+R4     11106.231    311 22834.463    35772.063    24013.140
203  TIM2+F+R5     11097.658    313 22821.316    37941.623    24007.573
204  TIM2+F+R6     11092.669    315 22815.339    40913.521    24009.176
209  TIM3e         13644.477    302 27892.954    35518.454    29037.522
210  TIM3e+I       11919.719    303 24445.438    32455.177    25593.796
211  TIM3e+G4      11263.301    303 23132.602    31142.341    24280.960
212  TIM3e+I+G4    11140.243    304 22888.487    31317.577    24040.634
213  TIM3e+R2      11457.862    304 23523.723    31952.814    24675.871
214  TIM3e+R3      11204.954    306 23021.907    32416.107    24181.635
215  TIM3e+R4      11136.523    308 22889.046    33463.713    24056.354
216  TIM3e+R5      11132.550    310 22885.100    34936.350    24059.987
222  TIM3+F        13780.112    305 28170.224    37058.796    29326.162
223  TIM3+F+I      11909.513    306 24431.026    33825.226    25590.754
224  TIM3+F+G4     11270.482    306 23152.965    32547.165    24312.692
225  TIM3+F+I+G4   11138.600    307 22891.201    32844.464    24054.719
226  TIM3+F+R2     11453.160    307 23520.321    33473.584    24683.839
227  TIM3+F+R3     11197.486    309 23012.972    34282.384    24184.070
228  TIM3+F+R4     11134.999    311 22891.998    35829.598    24070.675
229  TIM3+F+R5     11127.119    313 22880.238    38000.546    24066.496
230  TIM3+F+R6     11123.271    315 22876.541    40974.723    24070.379
235  TVMe          13701.428    303 28008.856    36018.595    29157.213
236  TVMe+I        11912.404    304 24432.809    32861.900    25584.957
237  TVMe+G4       11249.785    304 23107.571    31536.662    24259.719
238  TVMe+I+G4     11125.473    305 22860.946    31749.517    24016.884
239  TVMe+R2       11439.222    305 23488.444    32377.015    24644.382
240  TVMe+R3       11188.551    307 22991.103    32944.366    24154.621
241  TVMe+R4       11122.299    309 22862.598    34132.009    24033.695
242  TVMe+R5       11117.812    311 22857.624    35795.224    24036.301
248  TVM+F         13864.341    306 28340.682    37734.882    29500.410
249  TVM+F+I       11913.331    307 24440.662    34393.925    25604.180
250  TVM+F+G4      11249.086    307 23112.171    33065.434    24275.689
251  TVM+F+I+G4    11112.094    308 22840.189    33414.855    24007.496
252  TVM+F+R2      11433.603    308 23483.206    34057.873    24650.514
253  TVM+F+R3      11173.098    310 22966.197    35017.447    24141.085
254  TVM+F+R4      11113.042    312 22850.084    36800.941    24032.551
255  TVM+F+R5      11102.547    314 22833.094    39318.094    24023.141
256  TVM+F+R6      11098.304    316 22828.609    42863.009    24026.236
261  SYM           13641.639    304 27891.277    36320.368    29043.425
262  SYM+I         11909.737    305 24429.475    33318.046    25585.413
263  SYM+G4        11247.500    305 23104.999    31993.570    24260.937
264  SYM+I+G4      11124.667    306 22861.334    32255.534    24021.062
265  SYM+R2        11434.606    306 23481.212    32875.412    24640.940
266  SYM+R3        11188.369    308 22992.739    33567.405    24160.047
267  SYM+R4        11121.284    310 22862.568    34913.818    24037.456
268  SYM+R5        11116.242    312 22856.484    36807.341    24038.951
274  GTR+F         13773.461    307 28160.923    38114.186    29324.440
275  GTR+F+I       11890.831    308 24397.662    34972.329    25564.970
276  GTR+F+G4      11241.015    308 23098.030    33672.697    24265.338
277  GTR+F+I+G4    11104.880    309 22827.760    34097.172    23998.858
278  GTR+F+R2      11414.766    309 23447.532    34716.944    24618.630
279  GTR+F+R3      11165.214    311 22952.427    35890.027    24131.105
280  GTR+F+R4      11102.404    313 22830.807    37951.115    24017.065
281  GTR+F+R5      11095.292    315 22820.584    40918.766    24014.422
282  GTR+F+R6      11089.477    317 22812.954    45214.287    24014.371
283  GTR+F+R7      11089.312    319 22816.623    51982.338    24025.621
Akaike Information Criterion:           GTR+F+R6
Corrected Akaike Information Criterion: K2P+G4
Bayesian Information Criterion:         TIM2+F+I+G4
Best-fit model: TIM2+F+I+G4 chosen according to BIC

All model information printed to H3-aligned-mafft.fasta.model.gz
CPU time for ModelFinder: 84.980 seconds (0h:1m:24s)
Wall-clock time for ModelFinder: 85.155 seconds (0h:1m:25s)

NOTE: 4 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
Thoroughly optimizing +I+G parameters from 10 start values...
Init pinv, alpha: 0.000, 1.000 / Estimate: 0.000, 0.325 / LogL: -11242.309
Init pinv, alpha: 0.058, 1.000 / Estimate: 0.506, 0.805 / LogL: -11108.441
Init pinv, alpha: 0.116, 1.000 / Estimate: 0.506, 0.806 / LogL: -11108.443
Init pinv, alpha: 0.173, 1.000 / Estimate: 0.506, 0.805 / LogL: -11108.443
Init pinv, alpha: 0.231, 1.000 / Estimate: 0.506, 0.804 / LogL: -11108.444
Init pinv, alpha: 0.289, 1.000 / Estimate: 0.506, 0.805 / LogL: -11108.447
Init pinv, alpha: 0.347, 1.000 / Estimate: 0.506, 0.805 / LogL: -11108.447
Init pinv, alpha: 0.404, 1.000 / Estimate: 0.506, 0.804 / LogL: -11108.456
Init pinv, alpha: 0.462, 1.000 / Estimate: 0.506, 0.806 / LogL: -11108.429
Init pinv, alpha: 0.520, 1.000 / Estimate: 0.506, 0.806 / LogL: -11108.411
Optimal pinv,alpha: 0.506, 0.806 / LogL: -11108.411

Parameters optimization took 2.192 sec
Computing ML distances based on estimated model parameters... 0.128 sec
Computing BIONJ tree...
0.008 seconds
Log-likelihood of BIONJ tree: -11143.827
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 0.879 second
Computing log-likelihood of 98 initial trees ... 1.523 seconds
Current best score: -11087.389

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -11032.360
Iteration 10 / LogL: -11072.458 / Time: 0h:0m:6s
Iteration 20 / LogL: -11088.534 / Time: 0h:0m:8s
Finish initializing candidate tree set (20)
Current best tree score: -11032.360 / CPU time: 5.711
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 21: -11030.962
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 23: -11030.507
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 29: -11017.965
Iteration 30 / LogL: -11030.453 / Time: 0h:0m:10s (0h:0m:34s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 33: -11014.712
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 35: -11013.329
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 39: -11011.781
Iteration 40 / LogL: -11018.766 / Time: 0h:0m:12s (0h:0m:31s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 46: -11009.551
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 48: -11008.896
Iteration 50 / LogL: -11010.369 / Time: 0h:0m:14s (0h:0m:29s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 54: -11008.266
Iteration 60 / LogL: -11009.660 / Time: 0h:0m:16s (0h:0m:26s left)
Iteration 70 / LogL: -11014.711 / Time: 0h:0m:18s (0h:0m:23s left)
Iteration 80 / LogL: -11010.339 / Time: 0h:0m:20s (0h:0m:19s left)
Iteration 90 / LogL: -11009.090 / Time: 0h:0m:22s (0h:0m:16s left)
Iteration 100 / LogL: -11008.970 / Time: 0h:0m:25s (0h:0m:13s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 104: -11007.114
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 105: -11006.738
Iteration 110 / LogL: -11013.964 / Time: 0h:0m:27s (0h:0m:23s left)
Iteration 120 / LogL: -11010.291 / Time: 0h:0m:29s (0h:0m:20s left)
Iteration 130 / LogL: -11012.043 / Time: 0h:0m:31s (0h:0m:18s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 140: -11005.533
Iteration 140 / LogL: -11005.533 / Time: 0h:0m:33s (0h:0m:23s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 145: -11004.912
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 148: -11003.685
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 150: -11003.081
Iteration 150 / LogL: -11003.081 / Time: 0h:0m:35s (0h:0m:23s left)
Iteration 160 / LogL: -11023.829 / Time: 0h:0m:37s (0h:0m:21s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 169: -11002.568
Iteration 170 / LogL: -11022.463 / Time: 0h:0m:39s (0h:0m:23s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 176: -10999.070
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 179: -10996.866
Iteration 180 / LogL: -11001.482 / Time: 0h:0m:41s (0h:0m:22s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 181: -10996.557
Iteration 190 / LogL: -10997.975 / Time: 0h:0m:43s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 196: -10996.400
Iteration 200 / LogL: -11002.451 / Time: 0h:0m:45s (0h:0m:21s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 204: -10995.436
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 210: -10995.211
Iteration 210 / LogL: -10995.211 / Time: 0h:0m:47s (0h:0m:22s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 213: -10995.135
Iteration 220 / LogL: -10996.565 / Time: 0h:0m:49s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 224: -10995.003
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 229: -10994.851
Iteration 230 / LogL: -11001.108 / Time: 0h:0m:51s (0h:0m:22s left)
Iteration 240 / LogL: -10998.353 / Time: 0h:0m:53s (0h:0m:19s left)
Iteration 250 / LogL: -10996.019 / Time: 0h:0m:55s (0h:0m:17s left)
Iteration 260 / LogL: -11013.925 / Time: 0h:0m:57s (0h:0m:15s left)
Iteration 270 / LogL: -10995.337 / Time: 0h:0m:59s (0h:0m:12s left)
Iteration 280 / LogL: -10997.389 / Time: 0h:1m:1s (0h:0m:10s left)
Iteration 290 / LogL: -10996.579 / Time: 0h:1m:3s (0h:0m:8s left)
Iteration 300 / LogL: -10995.908 / Time: 0h:1m:5s (0h:0m:6s left)
Iteration 310 / LogL: -11011.326 / Time: 0h:1m:7s (0h:0m:4s left)
Iteration 320 / LogL: -11022.668 / Time: 0h:1m:8s (0h:0m:1s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 328: -10994.785
Iteration 330 / LogL: -10995.309 / Time: 0h:1m:10s (0h:0m:21s left)
Iteration 340 / LogL: -10996.640 / Time: 0h:1m:12s (0h:0m:18s left)
Iteration 350 / LogL: -10998.205 / Time: 0h:1m:14s (0h:0m:16s left)
Iteration 360 / LogL: -10996.236 / Time: 0h:1m:17s (0h:0m:14s left)
Iteration 370 / LogL: -10994.895 / Time: 0h:1m:19s (0h:0m:12s left)
BETTER TREE FOUND at iteration 380: -10994.781
Iteration 380 / LogL: -10994.781 / Time: 0h:1m:20s (0h:0m:21s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 388: -10992.833
Iteration 390 / LogL: -10995.901 / Time: 0h:1m:22s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 398: -10992.638
Iteration 400 / LogL: -10998.398 / Time: 0h:1m:24s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 409: -10991.606
Iteration 410 / LogL: -11023.612 / Time: 0h:1m:27s (0h:0m:21s left)
Iteration 420 / LogL: -11012.511 / Time: 0h:1m:29s (0h:0m:18s left)
Iteration 430 / LogL: -10993.103 / Time: 0h:1m:30s (0h:0m:16s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 431: -10989.215
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 440: -10988.604
Iteration 440 / LogL: -10988.604 / Time: 0h:1m:33s (0h:0m:21s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 447: -10987.654
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 449: -10986.530
Iteration 450 / LogL: -10988.146 / Time: 0h:1m:35s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 453: -10985.423
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 459: -10983.794
Iteration 460 / LogL: -10983.794 / Time: 0h:1m:37s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 465: -10983.715
Iteration 470 / LogL: -10988.491 / Time: 0h:1m:39s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 477: -10983.168
Iteration 480 / LogL: -10993.070 / Time: 0h:1m:41s (0h:0m:20s left)
Iteration 490 / LogL: -10983.893 / Time: 0h:1m:43s (0h:0m:18s left)
Iteration 500 / LogL: -10985.142 / Time: 0h:1m:45s (0h:0m:16s left)
Iteration 510 / LogL: -10984.856 / Time: 0h:1m:47s (0h:0m:14s left)
Iteration 520 / LogL: -10990.524 / Time: 0h:1m:49s (0h:0m:12s left)
Iteration 530 / LogL: -10986.703 / Time: 0h:1m:51s (0h:0m:9s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 534: -10983.045
Iteration 540 / LogL: -10983.201 / Time: 0h:1m:53s (0h:0m:19s left)
Iteration 550 / LogL: -11001.247 / Time: 0h:1m:55s (0h:0m:17s left)
Iteration 560 / LogL: -10983.235 / Time: 0h:1m:57s (0h:0m:15s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 563: -10982.555
Iteration 570 / LogL: -10986.051 / Time: 0h:1m:59s (0h:0m:19s left)
Iteration 580 / LogL: -10985.706 / Time: 0h:2m:1s (0h:0m:17s left)
Iteration 590 / LogL: -10983.253 / Time: 0h:2m:4s (0h:0m:15s left)
Iteration 600 / LogL: -10991.164 / Time: 0h:2m:6s (0h:0m:13s left)
Iteration 610 / LogL: -10984.057 / Time: 0h:2m:8s (0h:0m:11s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 611: -10982.127
Iteration 620 / LogL: -10983.487 / Time: 0h:2m:10s (0h:0m:19s left)
Iteration 630 / LogL: -11011.458 / Time: 0h:2m:12s (0h:0m:17s left)
Iteration 640 / LogL: -10992.934 / Time: 0h:2m:14s (0h:0m:14s left)
Iteration 650 / LogL: -10991.204 / Time: 0h:2m:16s (0h:0m:12s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 659: -10981.687
Iteration 660 / LogL: -10982.165 / Time: 0h:2m:18s (0h:0m:20s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 668: -10981.605
Iteration 670 / LogL: -10983.700 / Time: 0h:2m:20s (0h:0m:20s left)
Iteration 680 / LogL: -10984.015 / Time: 0h:2m:22s (0h:0m:18s left)
Iteration 690 / LogL: -10992.075 / Time: 0h:2m:24s (0h:0m:16s left)
Iteration 700 / LogL: -11007.260 / Time: 0h:2m:26s (0h:0m:14s left)
Iteration 710 / LogL: -10986.789 / Time: 0h:2m:28s (0h:0m:12s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 712: -10981.583
Iteration 720 / LogL: -11000.057 / Time: 0h:2m:31s (0h:0m:19s left)
Iteration 730 / LogL: -11014.929 / Time: 0h:2m:33s (0h:0m:17s left)
Iteration 740 / LogL: -10983.464 / Time: 0h:2m:35s (0h:0m:15s left)
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 748: -10981.553
Iteration 750 / LogL: -10987.302 / Time: 0h:2m:37s (0h:0m:20s left)
BETTER TREE FOUND at iteration 760: -10981.545
Iteration 760 / LogL: -10981.545 / Time: 0h:2m:39s (0h:0m:20s left)
Iteration 770 / LogL: -10981.667 / Time: 0h:2m:41s (0h:0m:18s left)
Iteration 780 / LogL: -10997.373 / Time: 0h:2m:43s (0h:0m:16s left)
Iteration 790 / LogL: -10984.533 / Time: 0h:2m:45s (0h:0m:14s left)
Iteration 800 / LogL: -10983.816 / Time: 0h:2m:47s (0h:0m:12s left)
Iteration 810 / LogL: -10981.664 / Time: 0h:2m:49s (0h:0m:10s left)
Iteration 820 / LogL: -10996.815 / Time: 0h:2m:51s (0h:0m:8s left)
Iteration 830 / LogL: -11011.796 / Time: 0h:2m:53s (0h:0m:6s left)
Iteration 840 / LogL: -10991.514 / Time: 0h:2m:55s (0h:0m:4s left)
Iteration 850 / LogL: -10987.351 / Time: 0h:2m:57s (0h:0m:2s left)
Iteration 860 / LogL: -10986.053 / Time: 0h:2m:59s (0h:0m:0s left)
TREE SEARCH COMPLETED AFTER 861 ITERATIONS / Time: 0h:2m:59s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -10981.545
Optimal log-likelihood: -10981.544
Rate parameters:  A-C: 1.89097  A-G: 4.18745  A-T: 1.89097  C-G: 1.00000  C-T: 5.86344  G-T: 1.00000
Base frequencies:  A: 0.230  C: 0.331  G: 0.274  T: 0.165
Proportion of invariable sites: 0.503
Gamma shape alpha: 0.766
Parameters optimization took 1 rounds (0.031 sec)
BEST SCORE FOUND : -10981.544
Total tree length: 11.914

Total number of iterations: 861
CPU time used for tree search: 177.189 sec (0h:2m:57s)
Wall-clock time used for tree search: 177.446 sec (0h:2m:57s)
Total CPU time used: 179.608 sec (0h:2m:59s)
Total wall-clock time used: 179.873 sec (0h:2m:59s)

Analysis results written to: 
  IQ-TREE report:                H3-aligned-mafft.fasta.iqtree
  Maximum-likelihood tree:       H3-aligned-mafft.fasta.treefile
  Likelihood distances:          H3-aligned-mafft.fasta.mldist
  Screen log file:               H3-aligned-mafft.fasta.log

Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/H3-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies/Gene_Trees

### Running ASTRAL
Caveat: Due to incredibly long run times for even a single gene in MrBayes, input for ASTRAL will rely on individual gene trees produced by IQ-Tree (Maximum Likelihood inference)
	Astral tends to be statistically inconsistent when gene trees are constructed by ML (strong impact of LBA), however with a focus on the interfamilial relationships in Assamiidae, I do not anticipate any significantly long branches

Creating Concatenated File for all 5 gene trees
1. <cd ~/Desktop/563-Final-Project/figures/IQTree-ML-Phylogenies
-moving to the subdirectory containing all 5 tree files from previous section
2. <cat *.treefile > input-genetrees-ML-cat.tre>
-cat will concatenate the five tree files into a single file named input-genetrees-ML-cat.tre
-the *.treefile specifies that all files containing the ".treefile" component in their name will be concatenated

ASTRAL
<java -jar ~/Desktop/Astral/astral.5.7.8.jar -i input-genetrees-ML-cat.tre -o output-speciestree-ASTRAL.tre 2> output-speciestree-ASTRAL.log>
-requires pathing to location of astral application, hence the ~/Desktop/Astral/astral5.7.8.jar
-i "filename" specifies the input file containing all your gene trees
-o "filename" creates the corresponding output file for the final tree
2> "filename.log" prints the Astral analysis info from the run to a permanent text file

output-speciestree-ASTRAL.tre moved to ./563-Final-Project/figures/ASTRAL/
input-genetrees-ML-cat.tre moved to ./563-Final-Project/figures/ASTRAL/
output-specestree-ASTRAL.log moved to ./563-Final-Project/analysis/ASTRAL-Outputs/

Visualizing the Output in FigTree
1. Once FigTree is loaded, select File -> Open, select output-speciestree-ASTRAL.tre
2. To make tree readable, slide expansion bar to increase separation between terminals
	-Under Appearance, increasing line weight
	-Under Tip Labels (should be checked), increase font size
3. Rooting
-select the branch with terminal "Troglosiro" 
	-a Cyphophthalmi (the most distantly related suborder of Opiliones)
-select reroot from the top menu
4. Support Values
-check node labels box; in display box, select label - nodes now labeled with likelihood scores

### Manuscript Figures
-treefiles visualized in FigTree before export to Adobe Illustrator for manual labelling of clades and taxa
		
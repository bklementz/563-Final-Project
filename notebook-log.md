# Description of Dataset
## Internal Phylogeny of Assamiidae
Assamiidae represents a highly diverse family of daddy longlegs (Arachnida: Opiliones) within the suborder Laniatores. They possess poorly resolved internal systematics based on historically bad species descriptions throughout the 1930's, making morphological classification difficult. 

### Included Taxa
Specimens were collected across Australasia and sampled from natural history collections including Harvard's Museum of Comparative Zoology and the Muséum d'histoire naturelle in Geneva. 59 sequenced specimens were included in subsequent analyses. 9 additional specimens were included from the dataset constructed by Giribet et al. 2010. 

Outgroup taxa represent sequences of two Beloniscidae and 11 Pyramidopidae. 

### Included Genes/Loci/Markers
DNA extraction, amplification, and sequencing yielded partial fragments of:
-two mitochondrial protein-coding genes: 16S rRNA and cytochrome c oxidase subunit I
-one nuclear protein-encoding gene: histone H3
-two nuclear ribosomal genes: 18S rRNA and 28S rRNA

# Quality Control
Dataset was provided cleaned and aligned from Prashant Sharma and Luciano Palmieri-Rocha. See the alignments under data/ in the folder "1_Sequence alignments".

Sequences were assembled and cleaned using Geneious 9.1.8.
Consensus sequences were submitted to BLAST against NCBI database to identify contaminations. 

# Alignment (Two Methods)
As sequences were provided already aligned, for the purposes of reproducing this step, alignment files were manually edited and all gaps removed from each sequence.

Edited sequences for 16S, 18S, 28S, COI, and H3 have been moved to the "Sequences-for-Alignment_Edited" in the ./data_clean subfolder.

## Mafft
### Installation
Mafft installed via this link: <https://mafft.cbrc.jp/alignment/software/macstandard.html> and selecting the mafft-7.490-signed.pkg
Upon download, I simply opened the .pkg file and followed the installer's instructions
### Description
Offers various multiple alignment strategies
1. FFT-NS-1, FFT-NS-2 (Progressive Methods)
-FFT-NS-1: first makes a rough distance matrix by counting number of shared 6-tuples between every sequence pair, then builds a tree, and finally aligns the sequences according to branching order
    -this is the simplest option and very fast.
-FFT-NS-2:
    -since the distance matrix is only approximated at the beginning, this method recomputes the guide tree from the FFT-NS-1 alignment and carries out a second progressive alignment from the recomputed tree

2. FFT-NS-i, NW-NS-i (Iterative Refinement Methods)
-Accuracy of progressive alignment improved by iterative refinement method. Implements simplified version of PRRN. 
-FFT-NS-i: initial alignment by FFT-NS-2 is subjected to iterative refinement, repeated until changes in WSP score cease or reaches 1,000 cycles.

3. L-INS-i, E-INS-i, G-INS-i (Iterative Refinement using WSP and Consistency Scores)
-Useful in difficult cases
-Use new objective function combining WSP score and COFFEE-like score to evaluate consistency between a myltiple alignment and pairwise alignment
-Pairwise alignments consist of three algorithms (N-W, local Smith-Waterman with affine gap costs, and local alignment with generalized affine gap costs)
-E-INS-i: useful in cases of conserved motifs in long unalignable regions
    -recommended if nature of sequences are not clear
    -Assumes that arrangement of conserved motifs shared by all sequences
-L-INS-i: useful to align set of sequences with sequences flanking one alignable domain
    -assumes that input sequences have only one alignable domain
G-INS-i: assumes entire region can be aligned and aligns them globally with N-W algorithm

### Limitations
-Accuracy: library extensions are not performed (as opposed to T-Coffee) because the authors belive that iterative refinement is more efficient than library extension
-Scalability: If two unrelated and long genomic DNA sequences are given, FFT-NS-2 model tries to make full-length alignment and uses a large CPU time. BLAST is a more suitable method in this case.
-Order of alignable blocks assumed to be conserved for all input sequences

### Running the Mafft Alignments
'cd Users/bklementz/Desktop/563-Final-Project/data_clean/Sequences-for-Alignment_Edited'
-moved into folder containing fasta files for alignment

### Alignment of 16S Sequences:
'mafft' - command to begin alignment
Input file? - 16S.fasta
Output file? - 16S-aligned-mafft.fasta
Output format? - 4. Fasta format, input order
Strategy - 5. L-INS-i (accurate)

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

 ### Alignment of 18S Sequences
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

 ### Alignment of 28S Sequences
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

 ### Alignment of COI Sequences
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

 ### Alignment of H3 Sequences
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

## MUSCLE
### Installation
Download link: https://drive5.com/muscle/downloads_v3.htm
Selected version for Mac OSX Intel i86 64 Bits "muscle3.8.31_i86darwin64.tar.gz"
Move to downloaded directory and untar with command <tar -zxvf muscle3.8.31_i86darwin64.tar.gz>
### Description
All info from Edgar 2004
-k mer distances computed for each pair of input sequences for constructing distance matrix
-Distance matrix used to produce binary tree
-progressive alignment following branching order of initial tree
-MUSCLE reestimates the tree using Kimura distance calculated from now aligned sequences (more accurate)
-these distance values (in a matrix), used to produce binary tree 2
-tree 2 is divided into two subtrees by deleting an edge, profile of multiple alignment in each subtree calculated
-If score is improved, new alignment is kept, otherwise discarded
-repeated until convergence or defined limit reached

### Performance
-MUSCLE outperformed T-Coffee, ClustalW, and Mafft in quality scores and CPU times; highest BAliBASE score
-better average accuracy and better speed than other software available at the time
-From Pervez et al 2014
  -Outperformed by ProbCons, SATe, MAFFT, and T-Coffee when evaluated based on the effect of indel size
  -similarly outperformed with increasing sequence size
  -performed the best with increasing deletion rates
  -T-Coffee, SATe also beat out MUSCLE when tested on BALiBASE benchmark alignments

### Running the MUSCLE Alignments
'cd Users/bklementz/Desktop/563-Final-Project/data_clean/Sequences-for-Alignment/'

### 16S Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in 16S.fasta -out 16S-aligned-muscle.fasta'
Name of Final Output File: "16S-aligned-muscle.fasta" in ./data_clean/Alignments/MUSCLE/

### 18S Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in 18S.fasta -out 18S-aligned-muscle.fasta'
Name of Final Output File: "18S-aligned-muscle.fasta" in ./data_clean/Alignments/MUSCLE/

### 28S Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in 28S.fasta -out 28S-aligned-muscle.fasta'
Name of Final Output File: "28S-aligned-muscle.fasta" in ./data_clean/Alignments/MUSCLE/

### COI Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in COI.fasta -out COI-aligned-muscle.fasta'
Name of Final Output File: "COI-aligned-muscle.fasta" in ./data_clean/Alignments/MUSCLE/

### H3 Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in H3.fasta -out H3-aligned-muscle.fasta'
Name of Final Output File: "H3-aligned-muscle.fasta" in ./data_clean/Alignments/MUSCLE/

Output Files for each alignment have been transferred to ~/563-Final-Project/data_clean/Alignments/MUSCLE/

# Distance and Parsimony Tree Construction
## Distance Methods
### Neighbor Joining
#### Description of Algorithm
-relies on an input of a distance matrix calculated based on the fraction of sites in which each pair of sequences differs (p-distance)
-the p-value differences are then converted into measures of evolutionary difference using different models
-Model relies on tree construction via Minimum Evolution
-minimum evolution will construct a tree based on minimizing the genetic differences between taxa (min. length)
-will construct additive, not ultrametric trees
#### Steps
1. Net divergence is calculated as the sum of distances for every end node in the dataset
2. a rate-corrected distance matrix is constructed from these values
3. a new node is defined that groups the two taxa that have the lowest rate-corrected distance value
4. Branch lengths from the new node to the daughter lineages are calculated. 
5. the branch lengths or genetic distance from the new node (grouping the two previous taxa) to the remaining taxa are then calculated
6. Steps are repeated until all taxa have been connected via a final node
#### Strengths
-will produce better fits to distance calculations when taxa experience unequal rates of evolution and follow non-molecular clock behavior
-highly scalable; remains fast to construct as number of taxa and sequence length increases
-statistically consistent under many models of evolution
  -will produce the correct tree topology as long as distance matrix is nearly additive
#### Limitations and Assumptions
-assumes unequal rates of evolution
-NJ examines only a tiny fraction of the total number of possible tree topologies
-may be biased by statistical error in distance measurement
  -estimates of distance (formulas or likelihood) produces standard errors that can be quite high unless sequences are very long
-loss of data; reduction of sequence data to single value distances
-can result in negative branch lengths (non-biological)

# Running the Distance Trees
All tree construction was done in RStudio
## Loading necessary packages
<library(ape)>
<library(adegenet)>
<library(phangorn)>

## Navigating working directory to folder containing alignment files
<getwd()> -determine which directory I'm currently in
<setwd('/Users/bklementz/Desktop/563-Final-Project/data/Alignments/Mafft')> -moving to my alignments directory

## Loading Sample Data (16S alignment)
<dna <- fasta2DNAbin(file="16S-aligned-mafft.fasta")>
-converts fasta alignment into a DNAbin object

## Computing Genetic Distances
-using TN93 model
    -allows for differing rates of transitions and transversion, heterogeneous base frequencies, and between-site varaition in substitution rates
<D <- dist.dna(dna, model="TN93")>

## Constructing NJ Tree
<tre <- nj(D)>
error: Error in nj(D) : missing values are not allowed in the distance matrix Consider using njs()
<tre <- njs(D)>

## Ladderizing the Tree
<tre <- ladderize(tre)>

## Plotting and Titling the Tree
<plot(tre, cex=.6)>
<title("Assamiidae-Phylo-16S-Mafft-Alignment)>

Following the same steps, NJ trees were constructed for each of the other four genes in the dataset.

Phylogeny could not be constructed using the 18S and 28S alignments. Produced an error message, "distance information insufficient to construct a tree, cannot calculate agglomeration criterion).

Phylogenies were exported as PDFs from RStudio to ~/563-Final-Project/figures/Distance-Phylogenies/


## Parsimony Methods
### Description
-relies on character-based data, either as one of 4 nucleotides or one of 20 amino acids in input sequences
-tree construction attempts to minimize the number of character state changes to explain the dataset
-with more character state changes, it implies a more complex hypothesis because homoplasy is an ad hoc hypothesis
### Steps
1. determine the number of character state changes needed to explain the data
2. search overall tree topologies and select the tree that minimizes the number of state changes
-Relies on the Fitch algorithm for equal costs of state changes, or the Sankoff algorithm for unequal costs.
### Fitch Algorithm
1. root tree in random place
2. calculate the set state of each internal node by..
  -forming the intersection of the two child state sets
  -if the intersection is non-zero, the set state for the node becomes the sum of the lengths of the child branches 
  -if the intersection is zero, the set state for the node becomes the sum of the lengths of the child branches plus one.
### Assumptions and Limitations
-relies on slow rates of evolution (not an assumption)
-independence among characters
-more computationally intensive than distance-based methods given it is an NP-hard problem
-can be heavily influenced by long branch attraction (statistical inconsistency)
-cannot account for convergent evolution
-not based on any model of evolution
-no guarantee it will construct the correct tree

# Running the Maximum Parsimony Trees
All tree construction was done in RStudio
## Loading necessary packages
<library(ape)>
<library(adegenet)>
<library(phangorn)>

## Input of sample data and conversion to phangorn object
<dna <- fasta2DNAbin(file="16S-aligned-mafft.fasta")>
<dna2 <- as.phyDat(dna)>

## Computing starting tree
-neighbor joining distance tree as initial topology for searching the tree space
-must compute parsimony score for this tree
<tre.ini <- nj(dist.dna(dna,model="raw"))>
<parsimony(tre.ini, dna2)>

## Searching for tree with maximum parsimony
<tre.pars <- optim.parsimony(tre.ini, dna2)>

## Plot tree
<plot(tre.pars, cex=0.6)>

MP trees for each of the other four genes in the dataset followed the same steps

Trees were exported from RStudio as pdfs and stored in the "~/563-Final-Project/figures/Parsimony-Phylogenies/" directory

### 16S
Parsimony score of initial tree: 4296
Parsimony score of final tree: 4091 after 43 nni operations

### 18S
Tree could not be constructed; error message following:
<tre.ini <- nj(dist.dna(dna,model="raw"))>
error: distance information insufficient to construct a tree, cannot calculate agglomeration criterion

### 28S
Tree could not be constructed; same error as above

### COI
Parsimony score of initial tree: 11477
Parsimony score of final tree: 11238 after 54 nni operations

### H3
Parsimony score of initial tree: 3201
Parsimony score of final tree: 2590 after 106 nni operations

# Maximum Likelihood Methods
## Chosen Method: IQ-Tree
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


### Limitations
-takes longer to run than RAxML when CPU time is variable
-IQ-Tree 1 only runs on unrooted trees
-primarily uses fast NNI, thought to produce less accurate trees than SPR
-stochastic NNI perturbations do not continue after 100 perturbations
-even with randomly sampled initial trees and perturbations, there is no guarantee the output phylogeny represents the global optimum.
-produces variation in log-likelihoods, can produce different tree topologies

## Running IQ-Tree
<cd ~/Desktop/563-Final-Project/data/Alignments/Mafft>
-navigating to directory containing my MSAs

<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s 16S-aligned-mafft.fasta>
-running 16S Mafft alignment 
-first specify location of iqtree application in downloaded directory
then -s <name of input file>
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/16S-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies

<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s 18S-aligned-mafft.fasta>
-running 18S Mafft alignment
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/18S-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies

<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s 28S-aligned-mafft.fasta>
-running 28S Mafft alignment
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/28S-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies

<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s COI-aligned-mafft.fasta>
-running COI Mafft alignment
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputss/COI-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies

<~/Desktop/iqtree-1.6.12-MacOSX/bin/iqtree -s H3-aligned-mafft.fasta>
-running H3 Mafft Alignment
Output files relocated to ~/563-Final-Project/analysis/IQTree-Outputs/H3-aligned-mafft-outputs/ and tree file transferred to ./figures/IQTree-ML-Phylogenies

# Bayesian Methods
## Chosen Method: MrBayes
### Installation
1. First install homebrew ("https://brew.sh/") using command </bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)">
2. From terminal window: 
	<brew tap brewsci/bio>
	<brew install mrbayes>
### Description of Algorithm
-Phylogenetic inferences based on posterior probability distribution of trees, calculated via Bayes theorem
	-the posterior distribution is proportional to the product of a likelihood function and a prior distribution (based on user's own understanding of the system)
-given the reliance on prior and likelihood functions/distributions, the posterior distribution is often intractable (overly complex)
-relies on Markov chain Monte Carlo (MCMC) to approximate posterior tree probabilities
-MCMC
	-new state of chain is proposed using stochastic process
	-acceptance probability is calculated
	-uniform random variable is selected
		-if less than the acceptance probability, new state is accepted and state of chain changes (will only accept changes that increase posterior probability)
		-if not, chain remains in same prior state
	-process is repeated thousands to millions of times
	-proportion of time a tree is visited during course of chain is an approximation of its posterior probability
	-uphill steps (over topology of posterior probability) are always accepted, slightly downhill steps are usually accepted, and drastic downhill steps are almost never accepted, thus chain tends to stay in regions of high posterior probability
-MCMC in Phylogenetic Context
	-start with a random tree and arbitrary values for branch lengths and model parameters
	-each generation consists of either proposing a new tree and accepting or rejecting the SPR, NNI, etc moves, or proposing and accepting/rejecting a new model parameter value
	-every k generations (user specified), the tree topology is saved (sample chain)
	-after n generations, it will summarize the sampled trees with histograms, means, credible intervals, and other analyses
	-once MCMC is completed, a list of all trees visited is produced and evaluate which tree was visited the most (where the algorithm spent the most time), suggesting higher posterior probability for that tree
-MrBayes can also run MCMCMC
	-runs a user defined number of chains, n-1 of these chains are heated
	-once all the chains have gone a step, swap attempted between two randomly selected chains
	-if accepted, chains switch states
	-heated chains more easily explore tree space (heating lowers peaks and fills valleys, smoother likelihood topology)
	-cold chain can leap valleys when swap occurs between cold chain (if stuck on local optimum) and heated chain exploring another peak
-Input data: aligned DNA or amino acid sequences in nexus format	
-MrBayes3
	-uses Metropolis-Hastings sampler and updates single parameters in each move
	-scales rate parameters so that branch lengths are based on the expected number of site changes

### Assumptions
-Likelihood function typically calculated under assumption that substitutions occur under time-homogeneous Poisson process
-"an important but commonly invoked constraint is the assumption of data homogeneity" - Ronquist and Huelsenbeck (2003)

### Strengths
-Advantages over other methods (Huelsenbeck and Ronquist 2001):
	-easy interpretation of results
	-incorporation of prior information
-offers standard MCMC algorithm, and Metropolis-coupled Markov chain Monte Carlo (MCMCMC)
-can be easily operated through command-line
-user can specify substitution model, prior distributions, details of MC analyses
-options to relax assumption of equal rates of evolution across sites (e.g., gamma-distributed)
-can infer ancestral states and accommodates uncertainty about tree/model parameters
-Markov chains based on Metropolis-Hastings algorithm more computationally efficient than ML boostrapping (Ronquist and Huelsenbeck 2003)
	-accelerates convergence of Markov chain
-datasets of more than 350 taxa require only moderate computational effort, tree spaces typically too large for ML bootstrapping
-MrBayes3 accommodates mixed models, data heterogeneity (AAs, nucleotides, morphology)
	-also offers single, doublet, and codon models for nucleotide data
	-AA data can be analyzed by fixed or variale rate matrices
-well suited to parallelization, can even yield linear speed-ups

### Limitations
-MrBayes3 - more parameters in the mixed model, visit each parameter more rarely, leads to slower mixing and a need to run the chains for longer before adequate sample of posterior distribution is reached
-does not support multithreading
-analyses often conducted with default priors; could yield biased/incorrect results (Nascimento et al. 2018)
-often slower than other phylogenetic methods, given traversal of tree space is hindered by keeping track of all the steps the chain has taken
-samples are often dependent due to small neighborhood of moves
	-requires thinning of the samples to approximate truly independent sampling

## Running MrBayes
### Converting MSAs from .fasta to .nxs
1. Opened file "COI-aligned-mafft.fasta" in ./data_clean/Alignments/Mafft/ with the alignment viewing software seaview
2. Selected File->Save As; in the new window, select Save As: NEXUS (.nxs)
3. Visual inspection of nexus file
	"BEGIN DATA;
		DIMENSIONS NTAX=120 NCHAR=657;
		FORMAT DATATYPE=DNA
		GAP=-
		;"
	and at the end...
		";
		END;"
	Matches characteristics of original .fasta file
4. "COI-aligned-mafft.nxs" saved to ./data_clean/Alignments/nex-format/
### Running MrBayes
1. Navigate to subdirectory containing nexus input file
	<cd Desktop/563-Final-Project/data_clean/Alignments/nexus-format>
2. Open MrBayes with command <mb>
3. <execute COI-aligned-mafft.nxs>
4. <lset nst=6 rates=invgamma> - selects GTR model with proportion of invariable sites and gamma distribution of rates across sites
5. <mcmc nruns=4 nchains=4 ngen=20000000 samplefreq=200>
		mcmc - starts the Markov chain Monte Carlo analysis
		nruns - sets how may independent analyses are started at once
		nchains - number of chains, default is 4 chains (1 cold and 3 hot)
		ngen - the number of generations analysis will run for (chose 20000000 based on advice from Prashant Sharma)
		samplefreq - how often the chain is sampled (=200 will result in 10000 samples over 20000000 generations)
6. Hitting enter after writing out the previous command begins the analysis
Initial Output:
MrBayes > execute COI-aligned-mafft.nxs

   Executing file "COI-aligned-mafft.nxs"
   UNIX line termination
   Longest line length = 61
   Parsing file
   Expecting NEXUS formatted file
   Reading data block
      Allocated taxon set
      Allocated matrix
      Defining new matrix with 120 taxa and 657 characters
      Data is Dna
      Gaps coded as -
      Taxon   1 -> Trojanella
      Taxon   2 -> Holoscotolemon
      Taxon   3 -> Peltonychia
      Taxon   4 -> Zuma
      Taxon   5 -> As009_Cam
      Taxon   6 -> As050_Cam
      Taxon   7 -> As018_Gab
      Taxon   8 -> As011_Gab
      Taxon   9 -> As056_Gab
      Taxon  10 -> As027_Gab
      Taxon  11 -> As084_Cam
      Taxon  12 -> As057_Gab
      Taxon  13 -> As058_Gab
      Taxon  14 -> As040_Gab
      Taxon  15 -> Scotolemon
      Taxon  16 -> Bishopella
      Taxon  17 -> Protolophus
      Taxon  18 -> Pantopsalis
      Taxon  19 -> Hesperonemastoma
      Taxon  20 -> Dendrolasma
      Taxon  21 -> Trogulus
      Taxon  22 -> Troglosiro
      Taxon  23 -> Synthetonychia
      Taxon  24 -> Conomma
      Taxon  25 -> Guasinia
      Taxon  26 -> Icaleptes_DNA104056
      Taxon  27 -> Ethobunus
      Taxon  28 -> Zalmoxis
      Taxon  29 -> Kimula
      Taxon  30 -> Stygnomma
      Taxon  31 -> Baculigerus
      Taxon  32 -> Stenostygnus_DNA104848
      Taxon  33 -> Stenostygnus_DNA104849
      Taxon  34 -> Stenostygnus_DNA104847
      Taxon  35 -> Stenostygnus_DNA104850
      Taxon  36 -> Fijicolana
      Taxon  37 -> Pellobunus
      Taxon  38 -> Metabiantes_DNA100704
      Taxon  39 -> Metabiantes_DNA100335
      Taxon  40 -> Biantidae_DNA105668
      Taxon  41 -> Metabiantes_DNA100703
      Taxon  42 -> Trionyxella
      Taxon  43 -> Paktongius
      Taxon  44 -> As120_Thai
      Taxon  45 -> As110_Laos
      Taxon  46 -> As121_Laos
      Taxon  47 -> As096_PNG
      Taxon  48 -> As095_PNG
      Taxon  49 -> As094_PNG
      Taxon  50 -> As102_NAus
      Taxon  51 -> As092_EAus
      Taxon  52 -> As104_EAus
      Taxon  53 -> As108_Laos
      Taxon  54 -> As081_Indo
      Taxon  55 -> As126_Laos
      Taxon  56 -> As080_Indo
      Taxon  57 -> As133_Indo
      Taxon  58 -> As127_Laos
      Taxon  59 -> As109_Thai
      Taxon  60 -> As087_WAus
      Taxon  61 -> Heterocranaus
      Taxon  62 -> Megapachylus
      Taxon  63 -> Goniosoma
      Taxon  64 -> Glysterus
      Taxon  65 -> Agoristenidae_DNA105839
      Taxon  66 -> Caenoncopus
      Taxon  67 -> Gnomulus
      Taxon  68 -> Palaeoncopus
      Taxon  69 -> Martensiellus
      Taxon  70 -> Sandokan
      Taxon  71 -> Epedanidae_DNA104066
      Taxon  72 -> Tithaeus
      Taxon  73 -> As106_Thai
      Taxon  74 -> As118_Thai
      Taxon  75 -> As115_Thai
      Taxon  76 -> As131_Mala
      Taxon  77 -> As114_Thai
      Taxon  78 -> As129_Laos
      Taxon  79 -> As122_Laos
      Taxon  80 -> As136_Viet
      Taxon  81 -> As123_Viet
      Taxon  82 -> As140_Viet
      Taxon  83 -> Epedanidae_DNA104062
      Taxon  84 -> Epedanidae_DNA104068
      Taxon  85 -> Zygopachylus
      Taxon  86 -> Cynortula
      Taxon  87 -> As026_Gab
      Taxon  88 -> Hoplobunus
      Taxon  89 -> Stygnopsidae_DNA103882
      Taxon  90 -> Karos
      Taxon  91 -> Stygnopsidae_DNA104855
      Taxon  92 -> Stygnopsidae_DNA104856
      Taxon  93 -> As020_Gab
      Taxon  94 -> As041_Gab
      Taxon  95 -> As030_Gab
      Taxon  96 -> Jarmilana
      Taxon  97 -> Arulla_DNA102666
      Taxon  98 -> As059_Gab
      Taxon  99 -> IC_DNA104070
      Taxon 100 -> Zalmoxida
      Taxon 101 -> As099_Dibunus
      Taxon 102 -> Op105_Dibunus
      Taxon 103 -> Op106_Toccolus
      Taxon 104 -> Op107_Nanepedanus_rufus
      Taxon 105 -> Op104_Dibuninae
      Taxon 106 -> Santinezia
      Taxon 107 -> As028_Gab
      Taxon 108 -> Lomanius_DNA104935
      Taxon 109 -> Santobius_DNA104931
      Taxon 110 -> As105_Phil
      Taxon 111 -> As083_Cam
      Taxon 112 -> Op049_Beloniscus
      Taxon 113 -> Bunofagea
      Taxon 114 -> As031_Gab
      Taxon 115 -> Fissiphallius
      Taxon 116 -> Larifuga
      Taxon 117 -> Triaenobunus
      Taxon 118 -> Triaenonychidae
      Taxon 119 -> Equitius
      Taxon 120 -> As119_Indo
      Successfully read matrix
      Setting default partition (does not divide up characters)
      Setting model defaults
      Seed (for generating default start values) = 1681318336
      Setting output file names to "COI-aligned-mafft.nxs.run<i>.<p|t>"
   Exiting data block
   Reached end of file

MrBayes > lset nst=6 rates=invgamma

   Setting Nst to 6
   Setting Rates to Invgamma
   Successfully set likelihood model parameters

MrBayes > mcmc nruns=4 nchains=4  ngen=2000000 samplefreq=200 

   WARNING: Reallocation of zero size attempted. This is probably a bug. Problems may follow.
   WARNING: Reallocation of zero size attempted. This is probably a bug. Problems may follow.
   Setting number of runs to 4
   Setting number of chains to 4
   Setting number of generations to 2000000
   Setting sample frequency to 200
   Setting chain output file names to "COI-aligned-mafft.nxs.run<i>.<p/t>"
   Running Markov chain
   MCMC stamp = 6096088174
   Seed = 284773038
   Swapseed = 1681318336
   Model settings:

      Data not partitioned --
         Datatype  = DNA
         Nucmodel  = 4by4
         Nst       = 6
                     Substitution rates, expressed as proportions
                     of the rate sum, have a Dirichlet prior
                     (1.00,1.00,1.00,1.00,1.00,1.00)
         Covarion  = No
         # States  = 4
                     State frequencies have a Dirichlet prior
                     (1.00,1.00,1.00,1.00)
         Rates     = Invgamma
                     The distribution is approximated using 4 categories.
                     Shape parameter is exponentially
                     distributed with parameter (1.00).
                     Proportion of invariable sites is uniformly dist-
                     ributed on the interval (0.00,1.00).

   Active parameters: 

      Parameters
      ---------------------
      Revmat              1
      Statefreq           2
      Shape               3
      Pinvar              4
      Ratemultiplier      5
      Topology            6
      Brlens              7
      ---------------------

      1 --  Parameter  = Revmat
            Type       = Rates of reversible rate matrix
            Prior      = Dirichlet(1.00,1.00,1.00,1.00,1.00,1.00)

      2 --  Parameter  = Pi
            Type       = Stationary state frequencies
            Prior      = Dirichlet

      3 --  Parameter  = Alpha
            Type       = Shape of scaled gamma distribution of site rates
            Prior      = Exponential(1.00)

      4 --  Parameter  = Pinvar
            Type       = Proportion of invariable sites
            Prior      = Uniform(0.00,1.00)

      5 --  Parameter  = Ratemultiplier
            Type       = Partition-specific rate multiplier
            Prior      = Fixed(1.0)

      6 --  Parameter  = Tau
            Type       = Topology
            Prior      = All topologies equally probable a priori
            Subparam.  = V

      7 --  Parameter  = V
            Type       = Branch lengths
            Prior      = Unconstrained:GammaDir(1.0,0.1000,1.0,1.0)


   Number of chains per MPI processor = 16

   The MCMC sampler will use the following moves:
      With prob.  Chain will use move
         0.93 %   Dirichlet(Revmat)
         0.93 %   Slider(Revmat)
         0.93 %   Dirichlet(Pi)
         0.93 %   Slider(Pi)
         1.85 %   Multiplier(Alpha)
         1.85 %   Slider(Pinvar)
         9.26 %   ExtSPR(Tau,V)
         9.26 %   ExtTBR(Tau,V)
         9.26 %   NNI(Tau,V)
         9.26 %   ParsSPR(Tau,V)
        37.04 %   Multiplier(V)
        12.96 %   Nodeslider(V)
         5.56 %   TLMultiplier(V)

   Division 1 has 507 unique site patterns
   Initializing conditional likelihoods

   Running benchmarks to automatically select fastest BEAGLE resource... A
B

   Using BEAGLE v4.0.0 (PRE-RELEASE) resource 0 for division 1:
      Rsrc Name : CPU (arm64)
      Impl Name : CPU-4State-Single
      Flags: PROCESSOR_CPU PRECISION_SINGLE COMPUTATION_SYNCH EIGEN_REAL
             SCALING_MANUAL SCALERS_RAW VECTOR_NONE THREADING_CPP      MODEL STATES: 4
   Initializing invariable-site conditional likelihoods

   Initial log likelihoods and log prior probs for run 1:
      Chain 1 -- -66586.032045 -- 165.615042
      Chain 2 -- -67015.101198 -- 165.615042
      Chain 3 -- -67508.013773 -- 165.615042
      Chain 4 -- -66091.164425 -- 165.615042

   Initial log likelihoods and log prior probs for run 2:
      Chain 1 -- -66862.308880 -- 165.615042
      Chain 2 -- -66712.047929 -- 165.615042
      Chain 3 -- -67185.087873 -- 165.615042
      Chain 4 -- -67550.534599 -- 165.615042

   Initial log likelihoods and log prior probs for run 3:
      Chain 1 -- -66575.987152 -- 165.615042
      Chain 2 -- -66720.220937 -- 165.615042
      Chain 3 -- -66250.851910 -- 165.615042
      Chain 4 -- -66807.916266 -- 165.615042

   Initial log likelihoods and log prior probs for run 4:
      Chain 1 -- -66501.983430 -- 165.615042
      Chain 2 -- -67399.220464 -- 165.615042
      Chain 3 -- -66459.410393 -- 165.615042
      Chain 4 -- -67161.511326 -- 165.615042


   Using a relative burnin of 25.0 % for diagnostics

   Chain results (2000000 generations requested):

Encountered several issues for troubleshooting...
1. The average standard deviation of split frequencies, after over 7 million generations, continued to increase, rather than converge towards 0 as expected
2. Feasibility of running analysis on my computer, after 7 million generations, estimated time to complete the total of 20 million generations was nearly 55 hours

Attempt 2
1. Using "COI-aligned-mafft.nxs" in ./data_clean/Alignments/nex-format
2. Creating an mbblock (series of commands to append to end of nexus document to run MrBayes)
	-in a new text file, named mbblock.txt
	"begin mrbayes;
	
		set autoclose=yes;
		lset nst=6 rates=invgamma;
		mcmcp Checkpoint=yes Checkfreq=2000000 ngen=20000000 printfreq=100000 samplefreq=5000 nchains=4 nruns=4 savebrlens=yes;
		mcmc;
		sumt nruns=4 Relburnin=yes burninfrac=0.20;
		
		end;"
		
	set autoclose=yes - when autoclose is set to yes, the program will not prompt me during execution
	mcmcp - mcmcp allows setting parameters before running the program
	checkpoint=yes - all current parameter values of chains will be printed to check-pointing file once reaching a Checkfreq generation, allows restarting of analysis from last checkpoint
	checkfreq=2000000 - will print current parameter values every 2000000 generations
	printfreq=100000 - will print info to screen every 100000 generations (recommended by coworkers to slightly reduce runtime)
	samplefreq=5000 - Markov chain sampled every 5000 generations, still results in 4000 samples (will increase independence)
	savebrlens - will save branch length info on trees
	sumt - will print summary tree at the end
	Relburnin=yes - specifies that a proportion of sampled values will be discarded as burnin
	burninfrac=0.20 - specifies that 20% of samples will be discarded
	
3. Append mbblock.txt to end of nexus file with <cat COI-aligned-mafft.nxs mbblock.txt > COI-mb.nxs
4. Run MrBayes with <mb COI-mb.nxs>
Initial Output:
Executing file "COI-mb.nxs"
   UNIX line termination
   Longest line length = 121
   Parsing file
   Expecting NEXUS formatted file
   Reading data block
      Allocated taxon set
      Allocated matrix
      Defining new matrix with 120 taxa and 657 characters
      Data is Dna
      Gaps coded as -
      Taxon   1 -> Trojanella
      Taxon   2 -> Holoscotolemon
      Taxon   3 -> Peltonychia
      Taxon   4 -> Zuma
      Taxon   5 -> As009_Cam
      Taxon   6 -> As050_Cam
      Taxon   7 -> As018_Gab
      Taxon   8 -> As011_Gab
      Taxon   9 -> As056_Gab
      Taxon  10 -> As027_Gab
      Taxon  11 -> As084_Cam
      Taxon  12 -> As057_Gab
      Taxon  13 -> As058_Gab
      Taxon  14 -> As040_Gab
      Taxon  15 -> Scotolemon
      Taxon  16 -> Bishopella
      Taxon  17 -> Protolophus
      Taxon  18 -> Pantopsalis
      Taxon  19 -> Hesperonemastoma
      Taxon  20 -> Dendrolasma
      Taxon  21 -> Trogulus
      Taxon  22 -> Troglosiro
      Taxon  23 -> Synthetonychia
      Taxon  24 -> Conomma
      Taxon  25 -> Guasinia
      Taxon  26 -> Icaleptes_DNA104056
      Taxon  27 -> Ethobunus
      Taxon  28 -> Zalmoxis
      Taxon  29 -> Kimula
      Taxon  30 -> Stygnomma
      Taxon  31 -> Baculigerus
      Taxon  32 -> Stenostygnus_DNA104848
      Taxon  33 -> Stenostygnus_DNA104849
      Taxon  34 -> Stenostygnus_DNA104847
      Taxon  35 -> Stenostygnus_DNA104850
      Taxon  36 -> Fijicolana
      Taxon  37 -> Pellobunus
      Taxon  38 -> Metabiantes_DNA100704
      Taxon  39 -> Metabiantes_DNA100335
      Taxon  40 -> Biantidae_DNA105668
      Taxon  41 -> Metabiantes_DNA100703
      Taxon  42 -> Trionyxella
      Taxon  43 -> Paktongius
      Taxon  44 -> As120_Thai
      Taxon  45 -> As110_Laos
      Taxon  46 -> As121_Laos
      Taxon  47 -> As096_PNG
      Taxon  48 -> As095_PNG
      Taxon  49 -> As094_PNG
      Taxon  50 -> As102_NAus
      Taxon  51 -> As092_EAus
      Taxon  52 -> As104_EAus
      Taxon  53 -> As108_Laos
      Taxon  54 -> As081_Indo
      Taxon  55 -> As126_Laos
      Taxon  56 -> As080_Indo
      Taxon  57 -> As133_Indo
      Taxon  58 -> As127_Laos
      Taxon  59 -> As109_Thai
      Taxon  60 -> As087_WAus
      Taxon  61 -> Heterocranaus
      Taxon  62 -> Megapachylus
      Taxon  63 -> Goniosoma
      Taxon  64 -> Glysterus
      Taxon  65 -> Agoristenidae_DNA105839
      Taxon  66 -> Caenoncopus
      Taxon  67 -> Gnomulus
      Taxon  68 -> Palaeoncopus
      Taxon  69 -> Martensiellus
      Taxon  70 -> Sandokan
      Taxon  71 -> Epedanidae_DNA104066
      Taxon  72 -> Tithaeus
      Taxon  73 -> As106_Thai
      Taxon  74 -> As118_Thai
      Taxon  75 -> As115_Thai
      Taxon  76 -> As131_Mala
      Taxon  77 -> As114_Thai
      Taxon  78 -> As129_Laos
      Taxon  79 -> As122_Laos
      Taxon  80 -> As136_Viet
      Taxon  81 -> As123_Viet
      Taxon  82 -> As140_Viet
      Taxon  83 -> Epedanidae_DNA104062
      Taxon  84 -> Epedanidae_DNA104068
      Taxon  85 -> Zygopachylus
      Taxon  86 -> Cynortula
      Taxon  87 -> As026_Gab
      Taxon  88 -> Hoplobunus
      Taxon  89 -> Stygnopsidae_DNA103882
      Taxon  90 -> Karos
      Taxon  91 -> Stygnopsidae_DNA104855
      Taxon  92 -> Stygnopsidae_DNA104856
      Taxon  93 -> As020_Gab
      Taxon  94 -> As041_Gab
      Taxon  95 -> As030_Gab
      Taxon  96 -> Jarmilana
      Taxon  97 -> Arulla_DNA102666
      Taxon  98 -> As059_Gab
      Taxon  99 -> IC_DNA104070
      Taxon 100 -> Zalmoxida
      Taxon 101 -> As099_Dibunus
      Taxon 102 -> Op105_Dibunus
      Taxon 103 -> Op106_Toccolus
      Taxon 104 -> Op107_Nanepedanus_rufus
      Taxon 105 -> Op104_Dibuninae
      Taxon 106 -> Santinezia
      Taxon 107 -> As028_Gab
      Taxon 108 -> Lomanius_DNA104935
      Taxon 109 -> Santobius_DNA104931
      Taxon 110 -> As105_Phil
      Taxon 111 -> As083_Cam
      Taxon 112 -> Op049_Beloniscus
      Taxon 113 -> Bunofagea
      Taxon 114 -> As031_Gab
      Taxon 115 -> Fissiphallius
      Taxon 116 -> Larifuga
      Taxon 117 -> Triaenobunus
      Taxon 118 -> Triaenonychidae
      Taxon 119 -> Equitius
      Taxon 120 -> As119_Indo
      Successfully read matrix
      Setting default partition (does not divide up characters)
      Setting model defaults
      Seed (for generating default start values) = 1681446703
      Setting output file names to "COI-mb.nxs.run<i>.<p|t>"
   Exiting data block
   Reading mrbayes block
      Setting autoclose to yes
      Setting Nst to 6
      Setting Rates to Invgamma
      Successfully set likelihood model parameters
      Setting check-pointing ('Checkpoint') to yes
      Setting check-pointing frequency to 2000000
      Setting number of generations to 20000000
      Setting print frequency to 100000
      Setting sample frequency to 5000
      Setting number of chains to 4
      WARNING: Reallocation of zero size attempted. This is probably a bug. Problems may follow.
      WARNING: Reallocation of zero size attempted. This is probably a bug. Problems may follow.
      Setting number of runs to 4
      Setting chain output file names to "COI-mb.nxs.run<i>.<p/t>"
      Successfully set chain parameters
      Running Markov chain
      MCMC stamp = 6171898991
      Seed = 705183611
      Swapseed = 1681446703
      Model settings:

         Data not partitioned --
            Datatype  = DNA
            Nucmodel  = 4by4
            Nst       = 6
                        Substitution rates, expressed as proportions
                        of the rate sum, have a Dirichlet prior
                        (1.00,1.00,1.00,1.00,1.00,1.00)
            Covarion  = No
            # States  = 4
                        State frequencies have a Dirichlet prior
                        (1.00,1.00,1.00,1.00)
            Rates     = Invgamma
                        The distribution is approximated using 4 categories.
                        Shape parameter is exponentially
                        distributed with parameter (1.00).
                        Proportion of invariable sites is uniformly dist-
                        ributed on the interval (0.00,1.00).

      Active parameters: 

         Parameters
         ---------------------
         Revmat              1
         Statefreq           2
         Shape               3
         Pinvar              4
         Ratemultiplier      5
         Topology            6
         Brlens              7
         ---------------------

         1 --  Parameter  = Revmat
               Type       = Rates of reversible rate matrix
               Prior      = Dirichlet(1.00,1.00,1.00,1.00,1.00,1.00)

         2 --  Parameter  = Pi
               Type       = Stationary state frequencies
               Prior      = Dirichlet

         3 --  Parameter  = Alpha
               Type       = Shape of scaled gamma distribution of site rates
               Prior      = Exponential(1.00)

         4 --  Parameter  = Pinvar
               Type       = Proportion of invariable sites
               Prior      = Uniform(0.00,1.00)

         5 --  Parameter  = Ratemultiplier
               Type       = Partition-specific rate multiplier
               Prior      = Fixed(1.0)

         6 --  Parameter  = Tau
               Type       = Topology
               Prior      = All topologies equally probable a priori
               Subparam.  = V

         7 --  Parameter  = V
               Type       = Branch lengths
               Prior      = Unconstrained:GammaDir(1.0,0.1000,1.0,1.0)


      Number of chains per MPI processor = 16

      The MCMC sampler will use the following moves:
         With prob.  Chain will use move
            0.93 %   Dirichlet(Revmat)
            0.93 %   Slider(Revmat)
            0.93 %   Dirichlet(Pi)
            0.93 %   Slider(Pi)
            1.85 %   Multiplier(Alpha)
            1.85 %   Slider(Pinvar)
            9.26 %   ExtSPR(Tau,V)
            9.26 %   ExtTBR(Tau,V)
            9.26 %   NNI(Tau,V)
            9.26 %   ParsSPR(Tau,V)
           37.04 %   Multiplier(V)
           12.96 %   Nodeslider(V)
            5.56 %   TLMultiplier(V)

      Division 1 has 507 unique site patterns
      Initializing conditional likelihoods

      Running benchmarks to automatically select fastest BEAGLE resource... A
B

      Using BEAGLE v4.0.0 (PRE-RELEASE) resource 0 for division 1:
         Rsrc Name : CPU (arm64)
         Impl Name : CPU-4State-Single
         Flags: PROCESSOR_CPU PRECISION_SINGLE COMPUTATION_SYNCH EIGEN_REAL
                SCALING_MANUAL SCALERS_RAW VECTOR_NONE THREADING_CPP         MODEL STATES: 4
      Initializing invariable-site conditional likelihoods

      Initial log likelihoods and log prior probs for run 1:
         Chain 1 -- -67239.229659 -- 165.615042
         Chain 2 -- -66682.355446 -- 165.615042
         Chain 3 -- -66184.689441 -- 165.615042
         Chain 4 -- -66927.925099 -- 165.615042

      Initial log likelihoods and log prior probs for run 2:
         Chain 1 -- -66849.459360 -- 165.615042
         Chain 2 -- -67047.518333 -- 165.615042
         Chain 3 -- -67141.276801 -- 165.615042
         Chain 4 -- -66901.636279 -- 165.615042

      Initial log likelihoods and log prior probs for run 3:
         Chain 1 -- -66856.402180 -- 165.615042
         Chain 2 -- -66661.747256 -- 165.615042
         Chain 3 -- -66268.938334 -- 165.615042
         Chain 4 -- -67128.740616 -- 165.615042

      Initial log likelihoods and log prior probs for run 4:
         Chain 1 -- -66057.261833 -- 165.615042
         Chain 2 -- -66773.780628 -- 165.615042
         Chain 3 -- -67177.095933 -- 165.615042
         Chain 4 -- -67331.081331 -- 165.615042


      Using a relative burnin of 25.0 % for diagnostics

		

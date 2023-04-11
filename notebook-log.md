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



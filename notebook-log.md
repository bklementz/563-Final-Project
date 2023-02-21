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

Edited sequences for 16S, 18S, 28S, COI, and H3 have been moved to the "Sequences-for-Alignment_Edited" in the "data" subfolder.

## Mafft
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
'cd Users/bklementz/Desktop/563-Final-Project/data/Sequences-for-Alignment_Edited'
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

Output Files for each alignment have been transferred to ~/563-Final-Project/data/Alignments/Mafft/

## MUSCLE
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
'cd Users/bklementz/Desktop/563-Final-Project/data/Sequences-for-Alignment/'

### 16S Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in 16S.fasta -out 16S-aligned-muscle.fasta'

### 18S Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in 18S.fasta -out 18S-aligned-muscle.fasta'

### 28S Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in 28S.fasta -out 28S-aligned-muscle.fasta'

### COI Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in COI.fasta -out COI-aligned-muscle.fasta'

### H3 Alignment
'~/Desktop/Muscle/muscle3.8.31_i86darwin64 -in H3.fasta -out H3-aligned-muscle.fasta'

Output Files for each alignment have been transferred to ~/563-Final-Project/data/Alignments/MUSCLE/


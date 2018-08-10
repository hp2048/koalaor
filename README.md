# koalaor
Code and data in this repository were used to identify OR genes in the koala genome.

# Requirements
1. [HMMER](http://hmmer.org/) for scanning genome sequence to find OR candidates (nhmmer) and removing false positives by scanning GPCR Class A (Rhodopsin like) HMM database (hmmscan). Please have `nhmmer` and `hmmscan` binaries in the PATH. 
2. [FASTA Package](https://github.com/wrpearson/fasta36)  for obtaining conceptual translations of OR genes. Please have `fasty36` in the PATH.
3. [MCL algorithm](https://micans.org/mcl/) to identify clusters of similar OR genes. This was used to select reference OR genes to construct HMM. This is not required to detect OR genes in a genome.

# Overview of the pipeline
```
nhmmer -o /dev/null --tblout nhmmer.output ORGenes.nt.hmm genome.fasta
perl extracthmmerhits.pl genome.fasta genome.orcandidates nhmmer.output
fasty36 -b 1 -z 11 -Q -d 1 -m "F10 fasty.output" genome.orcandidates ORGenes.aa.fa >/dev/null
perl getORORFs.pl fasty.output genome.orcandidates genome.output.OR class_a_rhodopsin_like.aa.hmm 1 speciedcode
```

Above commands were executed for each species examined on SGE cluster using `runORFinder.sh` script. I will tidy things up when a I get a chance. However, you should be able to use essential perl scripts and reference data provided here to "reproduce" our results.

# Methods in bit more detail

We downloaded 107,803 olfactory receptor CDS sequences in FASTA format from the NCBI Nucleotide Database by using '"olfactory receptor"[All Fields] AND "Chordata"[Organism] AND (biomol_mrna[PROP] AND refseq[filter])' as the search term on 22nd May 2015 along with corresponding “Feature Table". These sequences were filtered to obtain full-length OR sequences with uninterrupted open reading frame using extractCDS.pl script. Following filter steps were performed sequentially:
1. 4862 CDS without start codon at the beginning and stop codon at the end of the translation were removed. 
2. 26092 sequences without "olfactory receptor “ term in the FASTA header description were removed. 
3. 3991 sequences containing non-ACGT characters were removed.
4. 1672 sequences with non-amino acid characters were removed.
5. 4635 sequences shorter than 300aa length were removed.
8. 5679 sequences longer than 330aa were removed.

## Identification of clusters of similar OR sequences:
All-vsl-all BLASTP (version 2.2.30+, -evalue 1e-5) was performed using 60,872 remaining OR sequences that were retained after applying filters listed above. BLASTP alignments were processed to remove alignments to the self, alignments with percent identity <=50 and >=95 resulting in 27,858,657 pair-wise alignments of 60,842 OR amino acid sequences. Pair-wise sparse matrix of alignments between two OR amino acid sequences was created by using negative log10(e-value) as the distance measure. Format for the sparse matrix was `ORgene1 ORgene2 -log10(evalue)`. This sparse matrix was processed using MCL algorithm to find clusters of similar OR genes based on the e-value of their pair-wise alignments. We identified 102 clusters containing 2 to 4921 sequences within any given cluster. We selected two amino acid sequences (chosen randomly when more than 2 sequences present in a cluster) from each cluster and performed multiple sequence alignment by using  the Muscle software (v3.8.31, default parameters). Subsequently, Muscle was used to iteratively generate profile alignment of all 102 multiple sequence alignments. Amino-acid sequence alignment of 204 sequences were back translated to nucleotide alignments and HMM profile was build using the HMMER suit (version 3.1b1 May 2013). This HMM profile was subsequently used to scan reference genome sequences for identification of putative OR genes as described below.

## Identification of OR genes in reference genomes:
We chose genomes of 14 tetrapod species along with the Koala genome to identify and compare OR gene family repertoire. We used nhmmer program using the HMM profile for OR genes to identify candidate OR region in each of the genome assembly. Candidate regions of OR HMM hits along with 500nt flanking sequences were isolated and these sequences were aligned back to the 204 reference OR sequences using fasty tool from FASTASEARCH suit (version 36.8.8) to obtain conceptual translations. All fasty alignments were processed to identify open reading frames. Susequently, false positive sequences were removed by aligning all discovered candidate OR gene sequences to HMM database containing OR profile and other rhodopsin gene family profiles. If the best hit of a sequence was not the OR HMM then it was discarded as false positive and removed from further analysis.

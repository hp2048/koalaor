# koalaor
Code and data in this repository were used to identify OR genes in the koala genome.

# Requirements
1. [HMMER](http://hmmer.org/)
2. [FASTA Package for alignments](https://github.com/wrpearson/fasta36)

# Overview of the pipeline
```
nhmmer -o /dev/null --tblout nhmmer.output ORGenes.nt.hmm genome.fasta
perl extracthmmerhits.pl genome.fasta genome.orcandidates nhmmer.output
fasty36 -b 1 -z 11 -Q -d 1 -m "F10 fasty.output" genome.orcandidates ORGenes.aa.fa >/dev/null
perl getORORFs.pl fasty.output genome.orcandidates genome.output.OR class_a_rhodopsin_like.aa.hmm 1 speciedcode
```

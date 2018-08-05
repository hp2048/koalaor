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

Above commands were executed for each species examined on SGE cluster using `runORFinder.sh` script. I will tidy things up when a I get a chance. However, you should be able to use essential perl scripts and reference data provided here to "reproduce" our results.


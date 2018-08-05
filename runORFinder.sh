#!/bin/bash
set -Ex

####following parameters are passed as variables to qsub as follows
###qsub -v local_base=/state/partition1 species_nameid=Monodelphis_domestica.9315 -v nfs_base=/home/hardip/chordateORs -v threads=8 -v genome_file=GCF_000002295.2_MonDom5_genomic.fna.gz -v orrefaa=/home/hardip/olfactoryreceptors/reference/ORGenes.aa.fa -v orrefhmm=/home/hardip/olfactoryreceptors/reference/ORGenes.nt.hmm -v rhodospinhmm=/home/hardip/olfactoryreceptors/reference/rhodopsinfamily/class_a_rhodopsin_like.aa.hmm
echo `date`

function cleanup {
  rm -rf $local_base/$species_nameid/`basename $orrefaa`
  rm -rf $local_base/$species_nameid/`basename $orrefhmm`
  rm -rf $local_base/$species_nameid/`basename $genome_file .gz`
  rsync -a $local_base/$species_nameid/ $nfs_base/$species_nameid/
  rm -rf $local_base/$species_nameid/
  rm -rf $nfs_base/$species_nameid/`basename $genome_file .gz`.pipeline.running
  echo `date`
}
trap cleanup ERR INT TERM EXIT

execute_command () {
command="$1"
taskname="$2"
donefile="$3"
force_this=$4
if [ "$force_this" -eq 1 ] || [ ! -e $donefile ] || [ ! -s $donefile ] || [ "`tail -n1 $donefile | cut -f4 -d','`" != " EXIT_STATUS: 0" ]
then
force=1
if [ "$HOSTNAME" == "gduserv.anu.edu.au" ]
then
JOB_ID="localhost"
fi
eval /usr/bin/time --format='"COMMAND=%C\nELAPSED=%E, SECONDS=%e, CPU=%S, CPUPERCENT=%P, MAXRM=%M Kb, AVGRM=%t Kb, PAGEFAULTS=%F, RPAGEFAULTS=%R, SWAP=%W, WAIT=%w, FSI=%I, FSO=%O, SMI=%r, SMO=%s EXITSTATUS:%x"' -o $donefile -a -- $command
ret=$?
echo JOBID: $JOB_ID, TASKID: $SGE_TASK_ID, TASKNAME: $taskname, EXIT_STATUS: $ret >>$donefile
if [ "$ret" -ne 0 ]
then
echo ERROR: "$command"
echo WARN: $taskname failed with $ret exit code.
cleanup
exit $ret
fi
else
echo $taskname has finished with following details.
tail -n3 $donefile
fi
}


if [ -e $nfs_base/$species_nameid/`basename $genome_file .gz`.pipeline.done ]
then
	echo "pipeline has finished."
	exit 0
else

mkdir -p $local_base/$species_nameid
mkdir -p $nfs_base/$species_nameid

####copy data to local drive
command="rsync -a $nfs_base/$species_nameid/ $local_base/$species_nameid/"
execute_command "$command" copy2local $local_base/$species_nameid/`basename $genome_file .gz`.copy2local.done 1
command="rsync -a $orrefaa $local_base/$species_nameid/"
execute_command "$command" copyreforaa $local_base/$species_nameid/`basename $genome_file .gz`.copyreforaa.done 1
command="rsync -a $orrefhmm $local_base/$species_nameid/"
execute_command "$command" copyreforhmm $local_base/$species_nameid/`basename $genome_file .gz`.copyreforhmm.done 1

####unzip genome data
command="pigz -dkf -p $threads $local_base/$species_nameid/$genome_file"
execute_command "$command" unzip $local_base/$species_nameid/`basename $genome_file .gz`.unzip.done $force

####run NHMMER
command="nhmmer --cpu $threads -o /dev/null --tblout $local_base/$species_nameid/`basename $genome_file .gz`.nhmmer $local_base/$species_nameid/`basename $orrefhmm` $local_base/$species_nameid/`basename $genome_file .gz`"
execute_command "$command" nhmmer $local_base/$species_nameid/`basename $genome_file .gz`.nhmmer.done $force

####extract hits
command="perl /home/hardip/olfactoryreceptortools/extracthmmerhits.pl $local_base/$species_nameid/`basename $genome_file .gz` $local_base/$species_nameid/`basename $genome_file .gz`.orcandidates $local_base/$species_nameid/`basename $genome_file .gz`.nhmmer"
execute_command "$command" extracthits $local_base/$species_nameid/`basename $genome_file .gz`.extracthits.done $force

####run conceptual translation
if [ -s $local_base/$species_nameid/`basename $genome_file .gz`.orcandidates ]
then
command="fasty36 -b 1 -z 11 -Q -d 1 -T $threads -m \"F10 $local_base/$species_nameid/`basename $genome_file .gz`.fasty\" $local_base/$species_nameid/`basename $genome_file .gz`.orcandidates $local_base/$species_nameid/`basename $orrefaa` >/dev/null"
execute_command "${command}" fasty $local_base/$species_nameid/`basename $genome_file .gz`.fasty.done $force
else
echo No OR candidates found
exit 0
fi
###Annotate ORs
command="perl /home/hardip/olfactoryreceptortools/getORORFs.pl $local_base/$species_nameid/`basename $genome_file .gz`.fasty $local_base/$species_nameid/`basename $genome_file .gz`.orcandidates $local_base/$species_nameid/`basename $genome_file .gz` $rhodospinhmm $threads $speciescode"
execute_command "$command" getORF $local_base/$species_nameid/`basename $genome_file .gz`.getORF.done $force

touch $nfs_base/$species_nameid/`basename $genome_file .gz`.pipeline.done
command="rsync -a $local_base/$species_nameid/ $nfs_base/$species_nameid/"
execute_command "$command" copy2nfs $local_base/$species_nameid/`basename $genome_file .gz`.copy2nfs.done $force
rm -rf $nfs_base/$species_nameid/`basename $genome_file .gz`.pipeline.running
cleanup

fi
####touch $nfs_output_dir/pipeline.done

####http://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop
####above link has examples of parallelizing for loop in bash
####http://www.perlmonks.org/?node_id=791996
####link for Perl threads example


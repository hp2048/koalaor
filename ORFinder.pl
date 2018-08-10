#!/usr/bin/perl -l
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

##process arguments
my $config = configure(scalar @ARGV);
##check required binaries before progressing.
&checkbinaries;
##run nhmmer
runcommand("nhmmer --cpu $config->{'threads'} -o /dev/null --tblout $config->{'outputbase'}.nhmmer.out $config->{'orrefnthmm'} $config->{'reffasta'}");
##process nhmmer output
&processnhmmer;
##run fasty36 for conceptual translation
runcommand("fasty36 -b 1 -z 11 -Q -d 1 -m \"F10 $config->{'outputbase'}.OR.candidates.fasty.output\" $config->{'outputbase'}.OR.candidates.fa $config->{'orrefaafa'} >/dev/null");
##process fasty output
&processfasty;
##unzip rhodopsin_class_a hmm
runcommand("gunzip -c $config->{'rhodopsinhmm'}.gz > $config->{'rhodopsinhmm'}");
##run hmmscan to rule out false positives
runcommand("hmmscan --cpu $config->{'threads'} --tblout $config->{'outputbase'}.clean.vs.rhodopsinclassa.hmmout -o /dev/null $config->{'rhodopsinhmm'} $config->{'outputbase'}.firstpass.oraaclean.fa");
##process hmmscan output to generate final OR list
&processhmmscan;
##delete unwanted files
runcommand("rm -f $config->{'reffasta'} $config->{rhodopsinhmm}");


sub configure {
  my $args = shift;
	my $config = {};
	GetOptions(
	$config,
	'threads=i',
	'reffasta=s',
	'speciescode=s',
	'outputdir=s',
	) or usage(1);

	$config->{threads}      = 1 unless (defined($config->{threads}));
	$config->{orrefnthmm}   = abs_path(dirname(__FILE__)."/referencedata/ORGenes.nt.hmm");
	$config->{orrefaafa}    = abs_path(dirname(__FILE__)."/referencedata/ORGenes.aa.fa");
	$config->{geneticcode}  = getgeneticcode();
	$config->{rhodopsinhmm} = abs_path(dirname(__FILE__)."/referencedata/class_a_rhodopsin_family.aa.hmm");

  if ($config->{'outputdir'} && ! -e $config->{'outputdir'}) {
    runcommand("mkdir -p $config->{'outputdir'}");
    $config->{'outputdir'} = abs_path($config->{'outputdir'});
    $config->{'outputdir'} .= "/" unless ($config->{'outputdir'} =~ /\/$/);
  }
  else {
    print "\nERROR: Provide output directory path and ensure that it doesn't exist.";
    usage(1);
  }

  if ($config->{'reffasta'} && -e $config->{'reffasta'}) {
    $config->{'reffasta'} = abs_path($config->{'reffasta'});
    if ($config->{'reffasta'} =~ /\.gz$/){
      my $reffasta_base = basename($config->{'reffasta'}, ".gz");
      runcommand("gunzip -c $config->{'reffasta'} | head -n 5000000 > $config->{'outputdir'}$reffasta_base");
      $config->{'reffasta'} = "$config->{'outputdir'}$reffasta_base";
    }
  }
  else{
    print "\nERROR: Provide reference fasta sequence for detecting OR genes.";
    usage(1);
  }

  if ($config->{'speciescode'} && $config->{'speciescode'} !~ /\s/) {
    $config->{'outputbase'} = "$config->{'outputdir'}$config->{'speciescode'}";
  }
  else{
    print "\nERROR: Provide a species code ensuring no spaces in the code";
    usage(1);
  }

  return $config;
}


sub checkbinaries {
	##check for required binaries in the path
	runcommand("which nhmmer");
	runcommand("which gunzip");
	runcommand("which fasty36");
	runcommand("which hmmscan");

}

sub processhmmscan {
	my $hmm_hits = read_hmmscan("$config->{'outputbase'}.clean.vs.rhodopsinclassa.hmmout");
	my %hmm_hits = %$hmm_hits;

	my $aaseq = readFasta("$config->{'outputbase'}.firstpass.oraa.fa");
	my $ntseq = readFasta("$config->{'outputbase'}.firstpass.ornt.fa");


	open (NTOUT, ">$config->{'outputbase'}.OR.nt.final.fasta") or die $!;
	open (AAOUT, ">$config->{'outputbase'}.OR.aa.final.fasta") or die $!;
	open (GFF, ">$config->{'outputbase'}.OR.final.gff") or die $!;
	foreach my $tempor (keys %$ntseq){
		if (exists $hmm_hits{$tempor}){
			if ($hmm_hits{$tempor} eq "ChordateOR"){
				$tempor =~ /(\S+):(\d+):(\d+):([fr])/;
				my $chr = $1;
				my $start = $2;
				my $end = $3;
				my $strand = $4;
				my $status = "P";
				if (length($$aaseq{$tempor}{'sequence'})>=250 && substr($$aaseq{$tempor}{'sequence'},0,1) eq "M" && substr($$aaseq{$tempor}{'sequence'},-1) eq "*" && substr($$aaseq{$tempor}{'sequence'},0,-1) =~ /^[ACDEFGHIKLMNPQRSTVWY]+$/){
					$status = "F";
				}
				print NTOUT ">$config->{'speciescode'}:$chr:$start:$end:".(($strand eq "f") ? "1" : "-1").":$status\n$$ntseq{$tempor}{'sequence'}";
				print AAOUT ">$config->{'speciescode'}:$chr:$start:$end:".(($strand eq "f") ? "1" : "-1").":$status\n$$aaseq{$tempor}{'sequence'}";
				print GFF "$chr\tORFinder\tORGene\t$start\t$end\t.\t".(($strand eq "f") ? "+" : "-")."\t.\tID=$config->{'speciescode'}:$chr:$start:$end:".(($strand eq "f") ? "1" : "-1").":$status";
			}
			else{
				print "INFO: $tempor is a false positive for $config->{'speciescode'}";
				print "INFO: $hmm_hits{$tempor} is the HMM hit";
			}
		}
		else{
			print "INFO: No hits for $tempor in HMM" if (! exists $hmm_hits{$tempor});
		}
	}
}

sub read_hmmscan {
	my $hmmout    = shift;
	my %hmm_hits = ();
	open (HMM, "<$hmmout") or die $!;
	while (<HMM>){
		chomp $_;
		next if (substr($_,0,1) eq "#");
		my @a = split (/\s+/, $_);
		next if (exists $hmm_hits{$a[2]});
		$hmm_hits{$a[2]} = $a[0];
	}
	close HMM;
	return \%hmm_hits;
}

sub process_alignment {
  my $report = shift;
  my $qseq = shift;
  my %aln_info = ();
  my $hsp_counter = 0;
  my $hit_counter = 0;
  my $aln_cons = 0;
  my $current_subject_id = "";
  my $current_subject_description = "";
  foreach my $l (split(/\n/, $report)){
    ###skip up to the query ID
    next if ($l !~ /^>>>(\S+),.*(\d+) aa vs / && ! exists $aln_info{'query_id'});
    ###get query ID and length
    ###NOTE: may need to modify regular expression to accomodate sequence description in the future.
    if ($l =~ /^>>>(\S+),.*(\d+) aa vs /){
      $aln_info{'query_id'} = $1;
      $aln_info{'query_len'} = $2;
    }
    ###program info
    elsif ($l =~ /^; (pg\S+|mp\S+)?: (.*)/){
      $aln_info{$1} = $2;
    }
    ###get subject ID, description
    elsif ($aln_info{'query_id'} && $l =~ /^>>(\S+)\s*(.*)/){
      $hit_counter++;
      $hsp_counter=0;
      $current_subject_id = $1;
      $current_subject_description = $2;
      $aln_info{'Hits'}{$hit_counter}{'subject_id'} = $current_subject_id;
      $aln_info{'Hits'}{$hit_counter}{'subject_description'} = $current_subject_description;
    }
    ###get subject ID, description
    elsif ($aln_info{'query_id'} && $l =~ /^>--/){
      $hit_counter++;
      $hsp_counter=0;
      $aln_info{'Hits'}{$hit_counter}{'subject_id'} = $current_subject_id;
      $aln_info{'Hits'}{$hit_counter}{'subject_description'} = $current_subject_description;
    }
    ###hit info
    elsif ($l =~ /^; (f\S+|sw\S+)?: (.*)/){
      $aln_info{'Hits'}{$hit_counter}{$1} = $2;
    }
    ###mark hsps
    elsif ($hit_counter>0 && $l =~ /^>\S+/){
      $hsp_counter++;
      $aln_cons = 0;
    }
    ###alignment conservation sequence start
    elsif ($l =~ /^; al_cons:.*/){
      $aln_cons = 1;
    }
    ###get hsp query info
    elsif ($hsp_counter%2 != 0 && $l =~ /^; (sq\S+|al\S+)?: (.*)/){
      $aln_info{'Hits'}{$hit_counter}{'hsp'}{int($hsp_counter/2)+1}{"query_".$1} = $2;
    }
    ###get hsp subject info
    elsif ($hsp_counter%2 == 0 && $l =~ /^; (sq\S+|al\S+)?: (.*)/){
      $aln_info{'Hits'}{$hit_counter}{'hsp'}{int($hsp_counter/2)}{"subject_".$1} = $2;
    }
    ###query alignment sequence
    elsif ($hsp_counter%2 != 0 && $aln_cons == 0){
      $aln_info{'Hits'}{$hit_counter}{'hsp'}{int($hsp_counter/2)+1}{"query_alnseq"} .= $l;
    }
    ###subject alignment sequence
    elsif ($hsp_counter%2 == 0 && $aln_cons == 0){
      $aln_info{'Hits'}{$hit_counter}{'hsp'}{int($hsp_counter/2)}{"subject_alnseq"} .= $l;
    }
    ###alignment conservation sequence
    elsif ($aln_cons == 1){
      $aln_info{'Hits'}{$hit_counter}{'hsp'}{int($hsp_counter/2)}{"alncons_seq"} .= $l;
    }
    else{
      print "WARN: Should have not reached here.";
      print "WARN: ODD LINE: $l";
      print "WARN: QUERY= " . $aln_info{'query_id'} if ($aln_info{'query_id'});
      exit(1);
    }
  }
  my @orinfo = get_orinfo(\%aln_info, $qseq);
  return @orinfo;
}

sub get_orinfo {
  my $aln_info = shift;
  my $qseq = shift;
  my %aln_info = %$aln_info;
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}    = $aln_info{'Hits'}{1}{'hsp'}{1}{'query_alnseq'};
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}    =~ s/-//g;
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'}  = $aln_info{'Hits'}{1}{'hsp'}{1}{'query_al_start'};
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'}    = $aln_info{'Hits'}{1}{'hsp'}{1}{'query_al_stop'};
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} = $aln_info{'Hits'}{1}{'fy_frame'};
  ###add/remove 3 from the ORend because fasty output counts 1 AA short and hence 3nt short.
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} = ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f") ? $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} + 3 : $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} - 3;
  ###check last 10 positions for stop codon
  my $last10 = substr($aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}, -10);
  my $stop_position = index($last10, "*");
  my $end_length_difference = $aln_info{'Hits'}{1}{'hsp'}{1}{'subject_sq_len'} - $aln_info{'Hits'}{1}{'hsp'}{1}{'subject_al_stop'} + 1;
  if ($end_length_difference <= 10 && $stop_position >= 0 && $last10 !~ /[\/\\]/){
    $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'} = substr($aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}, 0, length($aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}) -  (10 - $stop_position) + 1);
    $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} = ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f") ? $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} - ((10 - $stop_position - 1) * 3) : $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} + ((10 - $stop_position - 1) * 3);
  }
  ###Find stop codon in the selected sequence
  else{
    my $extra_translation = "";
    if ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f"){
      for (my $i = $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'}; $i < length($$qseq{$aln_info{'query_id'}}{'sequence'}) ; $i+=3){
        my $codon = substr($$qseq{$aln_info{'query_id'}}{'sequence'}, $i, 3);
        $extra_translation .= "X" if ($codon =~ /N/);
        $extra_translation .= $config->{'geneticcode'}{$codon} if (exists $config->{'geneticcode'}{$codon});
        last if (substr($extra_translation, -1) eq "*");
      }
      if (substr($extra_translation, -1) eq "*"){
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} += length($extra_translation) * 3;
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'} .= $extra_translation;
      }
    }
    elsif ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "r"){
      for (my $i = $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} - 4; $i >= 0; $i-=3){
        my $codon = substr($$qseq{$aln_info{'query_id'}}{'sequence'}, $i, 3);
        $codon = reverse($codon);
        $codon =~ tr/[ACGT]/[TGCA]/;
        $extra_translation .= "X" if ($codon =~ /N/);
        $extra_translation .= $config->{'geneticcode'}{$codon} if (exists $config->{'geneticcode'}{$codon});
        last if (substr($extra_translation, -1) eq "*");
      }
      if (substr($extra_translation, -1) eq "*"){
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} -= length($extra_translation) * 3;
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'} .= $extra_translation;
      }
    }
  }

  ###check first 10 positions for start codon
  my $first10 = substr($aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}, 0, 15);
  my $start_position = index($first10, "M");
  my $start_length_difference = $aln_info{'Hits'}{1}{'hsp'}{1}{'subject_al_start'};

  if ($start_length_difference <= 10 && $start_position >= 0 && $first10 !~ /\/|\\/){
    $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'} = substr($aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'}, $start_position);
    $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} = ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f") ? $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} + ($start_position * 3) : $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} - ($start_position * 3);
  }
  ###Find start codon in the selected sequence
  else{
    my $extra_translation = "";
    if ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f"){
      for (my $i = $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} - 4; $i >= 0; $i-=3){
        my $codon = substr($$qseq{$aln_info{'query_id'}}{'sequence'}, $i, 3);
        $extra_translation = "X" . $extra_translation if ($codon =~ /N/);
        $extra_translation = $config->{'geneticcode'}{$codon} ."".$extra_translation if (exists $config->{'geneticcode'}{$codon});
        last if (substr($extra_translation, 0, 1) eq "*" || substr($extra_translation, 0, 1) eq "M");
      }
      if (substr($extra_translation, 0, 1) eq "M"){
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} -= length($extra_translation) * 3;
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'} = $extra_translation ."".$aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'};
        #print $$qseq{$aln_info{'query_id'}}{'sequence'};
        #die Dumper \%aln_info;
      }
    }
    elsif ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "r"){
      for (my $i = $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'}; $i < length($$qseq{$aln_info{'query_id'}}{'sequence'}); $i+=3){
        my $codon = substr($$qseq{$aln_info{'query_id'}}{'sequence'}, $i, 3);
        $codon = reverse($codon);
        $codon =~ tr/[ACGT]/[TGCA]/;
        $extra_translation = "X" . $extra_translation if ($codon =~ /N/);
        $extra_translation = $config->{'geneticcode'}{$codon} ."".$extra_translation if (exists $config->{'geneticcode'}{$codon});
        last if (substr($extra_translation, 0, 1) eq "*" || substr($extra_translation, 0, 1) eq "M");
      }
      if (substr($extra_translation, 0, 1) eq "M"){
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} += length($extra_translation) * 3;
        $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'} = $extra_translation ."".$aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'};
      }
    }
  }

  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORCDS'} = ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f") ? substr($$qseq{$aln_info{'query_id'}}{'sequence'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} - 1, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} - $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} + 1) :  substr($$qseq{$aln_info{'query_id'}}{'sequence'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} - 1, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} - $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} + 1) ;
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORCDS'} = reverse($aln_info{'Hits'}{1}{'hsp'}{1}{'ORCDS'}) if ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "r");
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORCDS'} =~ tr/[ACGT]/[TGCA]/ if ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "r");
  $aln_info{'query_id'} =~ /(\S+):(\d+)-(\d+)/;
  my $gname  = $1;
  my $gstart = $2;
  my $gend   = $3;
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigid'}    = $gname;
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigstart'} = ($aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'} eq "f") ? $gstart + $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstart'} - 1 : $gstart + $aln_info{'Hits'}{1}{'hsp'}{1}{'ORend'} - 1;
  $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigend'}   = $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigstart'} + length($aln_info{'Hits'}{1}{'hsp'}{1}{'ORCDS'}) - 1;
  die Dumper \%aln_info if ($aln_info{'query_id'} eq "ABVD02257015.1:3233-5802");
  #die Dumper \%aln_info if ($aln_info{'query_id'} eq "KN195678.1:245790-247744")
  return $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigid'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORstrand'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigstart'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORcontigend'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORCDS'}, $aln_info{'Hits'}{1}{'hsp'}{1}{'ORORF'};
}

sub processfasty {
	my %orgenes = ();
	my $queryseq = readFasta("$config->{'outputbase'}.OR.candidates.fa");
	my $aln_report = "";
	open (F, "<$config->{'outputbase'}.OR.candidates.fasty.output") or die $!;
	while (<F>){
		if ($_ =~ /\s*\d+>>>\S+/){
			if (length($aln_report)>0){
				my @orinfo = process_alignment($aln_report, $queryseq);
				$orgenes{$orinfo[0]}{$orinfo[2]}{$orinfo[3]}{$orinfo[1]}{$orinfo[4]}{$orinfo[5]}="";
			}
			$aln_report = $_;
		}
		elsif (length($aln_report) > 0) {
			$aln_report .= $_;
		}
	}
	close F;
	open (TNOUT, ">$config->{'outputbase'}.firstpass.ornt.fa");
	open (TAOUT, ">$config->{'outputbase'}.firstpass.oraa.fa");
	open (TACLEAN, ">$config->{'outputbase'}.firstpass.oraaclean.fa");
	foreach my $c (keys %orgenes){
		foreach my $s (keys %{$orgenes{$c}}){
			foreach my $e (keys %{$orgenes{$c}{$s}}){
				foreach my $t (keys %{$orgenes{$c}{$s}{$e}}){
					foreach my $n (keys %{$orgenes{$c}{$s}{$e}{$t}}){
						foreach my $a (keys %{$orgenes{$c}{$s}{$e}{$t}{$n}}){
							print TNOUT ">$c:$s:$e:$t\n$n";
							print TAOUT ">$c:$s:$e:$t\n$a";
							$a =~ s/[\\\/]//g;
							print TACLEAN ">$c:$s:$e:$t\n$a";
						}
					}
				}
			}
		}
	}
	close TNOUT;
	close TAOUT;
	close TACLEAN;
}

sub getgeneticcode {
	my %genetic_code  = (
	'TCA' => 'S','TCC' => 'S','AGC' => 'S','AGT' => 'S','TCG' => 'S','TCT' => 'S',
	'TTC' => 'F','TTT' => 'F',
	'TTA' => 'L','TTG' => 'L','CTC' => 'L','CTG' => 'L','CTT' => 'L','CTA' => 'L',
	'TAC' => 'Y','TAT' => 'Y',
	'TAA' => '*','TAG' => '*','TGA' => '*',
	'TGC' => 'C','TGT' => 'C',
	'TGG' => 'W',
	'CCA' => 'P','CCC' => 'P','CCG' => 'P','CCT' => 'P',
	'CAT' => 'H','CAC' => 'H',
	'CAA' => 'Q','CAG' => 'Q',
	'CGA' => 'R','CGC' => 'R','AGA' => 'R','AGG' => 'R','CGG' => 'R','CGT' => 'R',
	'ATA' => 'I','ATC' => 'I','ATT' => 'I',
	'ATG' => 'M',
	'ACA' => 'T','ACC' => 'T','ACG' => 'T','ACT' => 'T',
	'AAC' => 'N','AAT' => 'N',
	'AAA' => 'K','AAG' => 'K',
	'GTA' => 'V','GTC' => 'V','GTG' => 'V','GTT' => 'V',
	'GCA' => 'A','GCC' => 'A','GCG' => 'A','GCT' => 'A',
	'GAC' => 'D','GAT' => 'D',
	'GAA' => 'E','GAG' => 'E',
	'GGA' => 'G','GGC' => 'G','GGG' => 'G','GGT' => 'G',
	);
	return \%genetic_code;
}

sub processnhmmer {
	my $refseq = readFasta($config->{'reffasta'});
	open (OUT, ">$config->{'outputbase'}.OR.candidates.fa") or die $!;
	open (HMM, "<$config->{'outputbase'}.nhmmer.out") or die $!;
	while (<HMM>){
		chomp $_;
		next if (substr($_,0,1) eq "#");
		my @a = split (/\s+/, $_);
		if ($a[5] - $a[4] >= 500){
			my $start = ($a[11] eq "+") ? $a[8] - 501 : $a[9] - 501;
			my $end   = ($a[11] eq "+") ? $a[9] + 500 : $a[8] + 500;
			$start = 0 if ($start < 0);
			$end   = $a[10] if ($end > $a[10]);
			print OUT ">$a[0]:".($start+1)."-$end\n";
			print OUT uc(substr($$refseq{$a[0]}{'sequence'}, $start, $end - $start))."\n";
		}
	}
	close HMM;
	close OUT;
}

sub readFasta {
	my $file = shift;
	my %s = ();
	my $header = "";
	open (F, "<$file") or die $!;
	while (<F>){
		chomp $_;
		if ($_ =~ />(\S+)/){
			$header = $1;
		}
		else {
			$s{$header}{'sequence'} .= $_;
		}
	}
	close F;
	return \%s;
}

sub runcommand {
	my $command = shift;
	unless (system($command) == 0){
		print "ERROR: Follwing command failed";
		print $command;
		exit(1);
	}
}

sub usage {
	my $exit_code = shift;
  print <<USAGEMSG;

USAGE:

ORFinder.pl -reffasta mygenome.fa -speciescode ABCDE -outputdir outputdir

Options:
  -reffasta      reference fasta file for annotating OR genes. It can be gzipped.
  -speciescode   (Five) letter species code from UniProt Species Code
                 available at https://www.uniprot.org/docs/speclist.
                 NOTE: If your species code is not available then make
                 one using genus and species name. e.g. Xenopus tropicalis
                 code is XENTR. Don't include any spaces in the code.
	-outputdir     Path to output directory. It will be created by
	               the program. Program WILL NOT overwrite.
  -threads       number of threads to use (default = 1)
USAGEMSG
	exit($exit_code);
}

#!/usr/bin/perl
#

my $lacer_version = 0.425;

# version 0.425
# Replicate rows prior to svd to make sure we have more rows than columns
# This appears to be a new check in PDL somewhere between 2.007 and 2.018
#
# version 0.424
# change -dump to a flag instead of requiring an integer argument
# add percentage to # of recalibration bases on status display
# change default $MINCOV to 5 to prevent lots of coverage=1 bases taking up
# too many recalibration bases, particularly a problem when using -stopbases
#
# version 0.423
# clean up verbose vs. dump, default thread_debug is 0
# fixed bug with reporting reverse complement base when mapped to negative
# strand - context was ok, just the base at that position was revcomped when it
# shouldn't have been
#
# version 0.422
# add bed file input for regions
#
# version 0.421
# add some extra output for figures, quality control, etc.
#
# version 0.42
# fix up read group handling
# change defaults - no coverage limit
#
# version 0.41
# separate out GLOBALS - 0.5 experiment with consolidating made it slower
# don't use single anymore, use actual counts
# alter sort order to capture minor nonconsensus bases
#
# version 0.4
# start cleaning up a bit
# add a better progress indicator
#
# version 0.3
# thread enabled
# batch together pdl calls instead of doing one per base
#
# version 0.2
# change default range span to 0.8 instead of 0.5
# also make this a parameter
# add in context and cycle covariates
# switch histogram keying to enable this
#
use warnings;
use strict;
use Data::Dumper;
use threads;
use threads::shared;
use Thread::Queue;
use Bio::DB::Sam;
use PDL;
use PDL::NiceSlice;
use PDL::MatrixOps;
use PDL::Parallel::threads ("retrieve_pdls", "share_pdls", "free_pdls");
$PDL::BIGPDL = 1;
$PDL::BIGPDL = 1;	# so Perl doesn't warn
use Memory::Usage;
use Term::ProgressBar;
use Getopt::Long;
Getopt::Long::Configure("pass_through");
$SIG{INT} = \&sigint_handler;	# so we can stop a little early

# shared globals
my $ALLHIST = {}; share($ALLHIST);
my $VCFPOS = {}; share($VCFPOS);
my @RG_LIST = (); share(@RG_LIST);
my @THREADSTATUS = (); share(@THREADSTATUS); # started, data, goaway, interrupt
my @DUMPDATA = (); share(@DUMPDATA);
my ($CONSENSUS_TO_PRINT,
    $MAXCOV,
    $MINCOV,
    $MINMAPQ,
    $MINQ,
    $INCLUDEVCF,
    $MAXQUAL,
    $MINQUAL,
    $MINOR_FREQ,
    $DUMP,
    $ABORT,
    $CALIBRATION_BASES,
    $TOTAL_BASES,
    $USE_READGROUPS) :shared;
$CONSENSUS_TO_PRINT = 1;	# 0 means all
$MAXCOV = 0;	# max coverage at a given position - skip if greater
$MINCOV = 5;	# min coverage at a given position - skip if less
$MINMAPQ = 30;	# minimum mapping quality to consider
$MINQ = 6;		# minimum base quality to consider
$INCLUDEVCF = 0;	# 0 is skip vcf positions, 1 is only use vcf positions
$MAXQUAL = 0;
$MINQUAL = 0;
$MINOR_FREQ = 0.05;	# we'll try to print these bases - generalized "single"
$DUMP = 0;
$ABORT = 0;	# second control-C we really quit
$CALIBRATION_BASES = 0;	# these are for dynamic updates
$TOTAL_BASES = 0;
$USE_READGROUPS = 1;	# if it's 0 we put everything together

# other variables - shouldn't be needed in threads
my $original_command = join (" ", $0, @ARGV);
my $mu = Memory::Usage->new();
my $starttime = time();
my $bamfile = "";
my $ref_fasta = "";
my $windowsize = 10000;	# windows size of coordinates to look at
my $region = "";	# look only in this region
my $outfile = "-";
my $vcf = "";
my $vcfoffset = 1000;
my $svd_bin = 3000;
my $covariate_binsize = 1000;
my $min_span = 0.8;
my $do_covariates = 1;
my $gatk = 1;		# whether to output gatk recal file
my $rg = "";
my $verbose = 1;
my $num_threads = 4;
my $rgfield = "ID";
my $randomize_regions = 0;
my $stop_bases = 0;
my $thread_debug = 0;
GetOptions (
  'bam=s' => \$bamfile,
  'reference=s' => \$ref_fasta,
  'consensus=i' => \$CONSENSUS_TO_PRINT,
  'window=i' => \$windowsize,
  'svdbin=i' => \$svd_bin,
  'covariatebin=i' => \$covariate_binsize,
  'coverage=i' => \$MAXCOV,
  'mincoverage=i' => \$MINCOV,
  'span=f' => \$min_span,
  'mapq=i' => \$MINMAPQ,
  'minq=i' => \$MINQ,
  'region=s' => \$region,
  'vcf=s' => \$vcf,
  'includevcf=i' => \$INCLUDEVCF,
  'gatk!' => \$gatk,
  'outfile=s' => \$outfile,
  'randomize!' => \$randomize_regions,
  'stopbases=i' => \$stop_bases,
  'verbose=i' => \$verbose,
  'covariates!' => \$do_covariates,
  'threads=i' => \$num_threads,
  'minor=f' => \$MINOR_FREQ,
  'dump!' => \$DUMP,
  'readgroups!' => \$USE_READGROUPS,
  'rgfield=s' => \$rgfield
);

if (!-f $bamfile || !-f $ref_fasta) {
  print "Usage: $0 -bam <bam file> -reference <ref fasta> [ parameters ]\n";
  print <<__USAGE__;
Required arguments
  -bam <bam file>            : make sure it's *sorted* and *indexed*
                               make sure the index file is named the same as
                               your .bam file but with a .bai index
  -reference <ref fasta>     : the reference fasta file

General data parameters (mostly similar to GATK)
  -vcf <vcf file>            : vcf file marking known variants - generally not
                               needed (no default)
  -includevcf <0|1>          : 0 = skip vcf positions (like GATK)
                               1 = use only vcf positions
  -coverage <int>            : if coverage at a position is greater than this,
                               don't pick single/nonconsensus bases here
                               (0 = no limit) (default $MAXCOV)
  -mincoverage <int>         : if coverage at a position is less than this,
                               don't pick single/nonconsensus bases here
                               (0 = no minimum) (default $MINCOV)
  -minq <int>                : minimum base quality to use, lower ones ignored
                               (default $MINQ)
  -mapq <int>                : minimum mapping quality - discard read otherwise
                               (default $MINMAPQ)
  -region <chrom:first-last> : region to focus on for pileups and recalibration
                               can also specify a BED file (chrom first last)
                               BED format is whitespace delimited
                               (default all chromosomes, all positions)
  -covariates|nocovariates   : whether to do context and cycle covariates
                               (default covariates)
  -gatk|nogatk               : whether to output GATK-style recalibration table
                               (default gatk)
  -readgroups|noreadgroups   : whether to use readgroups - if noreadgroups,
                               all data in the bam file will be combined
                               (default readgroups)
  -rgfield                   : specifies whether to use read group ID (ID),
                               platform unit (PU), library (LB), or sample
                               name (SM) from the \@RG line (later versions of
                               GATK require PU) (default $rgfield)

General running parameters
  -outfile <filename>        : output filename; default is STDOUT
  -threads <int>             : # of threads to use (default $num_threads)
  -window <int>              : window size (work unit) for pileups ($windowsize)
  -randomize|norandomize     : randomize order of windows for pileups
                               (default norandomize)
  -stopbases <int>           : stop pileups and recalibrate when we reach this
                               many recalibration bases (0 = no limit)
                               (default $stop_bases)

Technical SVD/recalibration parameters
  -consensus <int>           : # of consensus bases to use
                               every time we see a nonconsensus base, we take
                               this many consensus bases from that same pos
                               (0 = use all) (default $CONSENSUS_TO_PRINT)
  -svdbin <int>              : # of bases to bin into quality histogram when
                               making matrix that we end up doing SVD on ($svd_bin)
  -covariatebin <int>        : # of bases to bin when doing covariate SVDs
                               covariates generally have less data so make this
                               smaller than svdbin ($covariate_binsize)
  -span <float>              : between 0 and 1, technical parameter for how much
                               of eigenvalue 1 range to require when predicting
                               error and correct quality dists (default $min_span)
  -minor <float>             : between 0 and 1, realistically less than 0.1
                               defines frequency of minor bases in a pileup
                               (default $MINOR_FREQ)
  -verbose <0|1|2>           : 0 = totally silent
                               1 = progress meter (default)
                               2 = noisy (memory usage, timing, etc on STDOUT)
  -dump|nodump               : dump all the data to STDOUT (warning it's a lot)
__USAGE__
  exit;
}

# set up variables and files and such
if (!length($outfile) || $outfile eq "-") {
  *OUT = *STDOUT;
} elsif (-f $outfile) {
  die "Output file $outfile exists, will not overwrite.  Aborting...\n";
} else {
  open(OUT, ">$outfile");
}
my $sam = Bio::DB::Sam->new(-fasta => $ref_fasta, -bam => $bamfile, -autoindex => 1, -expand_flags => 1);
my $bam = $sam->bam;
my $header = $bam->header();
my $bam_comments = $header->text;
my @seqs = $sam->seq_ids;
my $i = 0;
my $first = 0;
my $last = 0;
my $single = ();
my $alldata = {};
my $threaddata = [];
my $tempdata;
my $temphist;
my @h = ();
my @regions = ();
my @out;
my %context_code = (
  "AA" =>  1, "AC" =>  2, "AG" =>  3, "AT" =>  4,
  "CA" =>  5, "CC" =>  6, "CG" =>  7, "CT" =>  8,
  "GA" =>  9, "GC" => 10, "GG" => 11, "GT" => 12,
  "TA" => 13, "TC" => 14, "TG" => 15, "TT" => 16,
  ".A" => -1, ".C" => -2, ".G" => -3, ".T" => -4,
  "A." => -5, "C." => -6, "G." => -7, "T." => -8,
  ".." => -9
);
my %code_context = reverse %context_code;
my $align;
my @f;
my $j;
my $time;
my $end;
my $last_count;
my $matrix;
my $recalibrated;
my $scale;
my $range;
my $chrom;
my @data;
my $p;
my $ix;
my $start;
my $consensus_base;
my $svd_hist;
my $covariate;
my $covariate_hist;
my $covariate_recalibrated;
my $queue;
my $number_regions;
my @threadlist;
my $temppdl;
my $pending;
my $waiting;
my $wait_interval;
my $progress;
my @table;
my $max_bases;
my $rginfo;
my $left_fields;

$verbose = 2 if $verbose > 2;
if ($verbose == 2) {
  print STDOUT "# Start time: $starttime\n";
  print STDOUT "# Original command line: ", $original_command, "\n";
  print STDOUT "# windowsize (for pileups): $windowsize\n";
  print STDOUT "# SVD binsize (for quality histograms): $svd_bin\n";
  print STDOUT "# Maxcov $MAXCOV, Min mapping $MINMAPQ, Min base $MINQ\n";
  print STDOUT "# vcf mode $INCLUDEVCF\n";
  $mu->record("Start");
}

#
# figure out what regions we're looking at
#
if (length $region && -f $region) {
  open REG, $region;
  while ($j = <REG>) {
    next if $j =~ /^#/;
    next if $j =~ /^$/;
    chomp;
    ($chrom, $start, $end) = split /\s+/, $j;
    # bed file is 0-based so correct range, but only need start
    $start++;
    for ($i = $start; $i <= $end; $i += $windowsize) {
      $first = $i;
      $last = $i+$windowsize-1;
      $last = $end if $end < $last;
      push @regions, "$chrom:$first-$last";
    }
  }
} elsif (length $region) {
  ($chrom, $range) = split /:/, $region;
  ($start, $end) = split /-/, $range;
  if (!$start || !$end) {
    $start = 1;
    $end = $sam->length($chrom);
  }
  for ($i = $start; $i <= $end; $i += $windowsize) {
    $first = $i;
    $last = $i+$windowsize-1;
    $last = $end if $end < $last;
    push @regions, "$chrom:$first-$last";
  }
}
if (!scalar @regions) {
  foreach $chrom (@seqs) {
    $start = 1;
    $end = $sam->length($chrom);
    for ($i = $start; $i <= $end; $i += $windowsize) {
      $first = $i;
      $last = $i+$windowsize-1;
      $last = $end if $end < $last;
      push @regions, "$chrom:$first-$last";
    }
  }
}
if ($randomize_regions) {
  my @temp = @regions;
  @regions = ();
  my @index = ();
  foreach $i (0..$#temp) {
    $index[$i] = rand();
  }
  foreach $i (sort { $index[$a] <=> $index[$b] } (0..$#temp)) {
    push @regions, $temp[$i];
  }
}

#
# hope the read groups are all listed in the bam comments at the top
# use this to initialize the pdls for $alldata
#
if ($USE_READGROUPS) {
  undef $rginfo;
  @f = split /\n/, $bam_comments;
  foreach $i (@f) {
    if ($i =~ /^\@RG/) {
      $i =~ /ID:(\S+)/;
      $rg = $1;
      $rg = "NULL" if !defined $rg;
      push(@RG_LIST, $rg);
      foreach $j (qw(ID PL PU LB SM)) {
        if ($i =~ /$j:(\S+)/) {
          $rginfo->{$rg}->{$j} = $1;
        } else {
          $rginfo->{$rg}->{$j} = "NULL";
        }
      }
    }
  }
  if (!defined $rginfo) {
    @RG_LIST = ("NULL");
    foreach $j (qw(ID PL PU LB SM)) {
      $rginfo->{"NULL"}->{$j} = "NULL";
    }
  }
}
if (!scalar @RG_LIST) {
  @RG_LIST = ("NULL");
}
foreach $rg (@RG_LIST) {
  foreach $j (1..$num_threads) {
    $threaddata->[$j]->{$rg} = zeroes(6,1);
    share_pdls($j . "__" . $rg => $threaddata->[$j]->{$rg});
  }
  $alldata->{$rg} = null;
}

#
# generate pileups and collect info
#
if (-f $vcf) {
  open V, $vcf;
  while (<V>) {
    next if /^#/;
    chomp;
    @f = split /\t/, $_;
    if (!defined $VCFPOS->{$f[0]}) {
      my %anon :shared;
      $VCFPOS->{$f[0]} = \%anon;
    }
    $VCFPOS->{$f[0]}->{$f[1]} = 1;
  }
  close V;
}

$queue = Thread::Queue->new;
$queue->enqueue(@regions);
$number_regions = scalar @regions;

$num_threads = scalar @regions if scalar @regions < $num_threads;
foreach $i (1..$num_threads) {
  threads->create(\&worker, $queue, $i, $i . "__");
  $THREADSTATUS[$i] = "started";
}

sub worker {
  my $q = shift;
  my $thread = shift;
  my $prefix = shift;
  my $region = "";
  my $time;
  my $tempdata;
  my $temphist;
  my $rg;
  my $first;
  my $last;
  my $alldata = {};
  foreach $rg (@RG_LIST) {
    $alldata->{$rg} = retrieve_pdls($prefix . $rg);
  }
  my $map_consensus = sub {
    my ($seqid, $pos, $pileup, $sam) = @_;
    my @a;
    my $aln;
    my @nt;
    my @rg;
    my $cov = 0;
    my $max = 0;
    my @con = ();
    my $num = 0;
    my $tempq;
    my $temppos;
    my @index;
    my $left;
    my $left2;
    my $right;
    my $printed_consensus;
    my $printed_nonminor;
    my @context;
    my $i;
    my $n;
    my $base = { "G" => 0, "A" => 0, "T" => 0, "C" => 0 };
    my $have_minor = 0;

    # first pass - call consensus
    return if ($pos < $first || $pos > $last);
    if ($INCLUDEVCF) {
      return if !defined $VCFPOS->{$seqid}->{$pos};
    } else {
      return if defined $VCFPOS->{$seqid} && defined $VCFPOS->{$seqid}->{$pos};
    }
    # keep track of indices in original @$pileup
    # only want to run alignments once
    # also need to keep track if we skipped this element of @$pileup
    # so initialize to 0
    @a = (0) x scalar(@$pileup);
    # do reverse so we can splice
    foreach $i (reverse (0..$#{$pileup})) {
      $p = $pileup->[$i];
      next if !defined($p);
      # if we don't want this read then don't look at it again
      next if ($p->is_refskip || $p->indel);

      $a[$i] = $p->alignment;
      $aln = $a[$i];
      next if !$aln->qual || $aln->qual < $MINMAPQ;
      $tempq = unpack("x". ($p->qpos) . "C", $aln->_qscore);
      next if $tempq < $MINQ;
      $cov++;
      if ($USE_READGROUPS) {
        $rg[$i] = $aln->aux;
        if ($rg[$i] =~ /RG:Z:(\S+)/) {
          $rg[$i] = $1;
        } else {
          $rg[$i] = "NULL";
        }
      } else {
        $rg[$i] = "NULL";
      }
      $temphist->{$rg[$i]}->{0}->[$tempq]++;	# overall histogram
      if ($aln->strand > 0) {
        # note qpos returns 0-based coordinate
        # context should be reverse complemented already if negative strand
        # base returned should be positive strand regardless - to compare with ref
        $temppos = $p->qpos+1;
#        ($context[$i], $nt[$i]) = get_context($aln->qseq, $p->qpos, 1);
        ($context[$i], $nt[$i]) = get_2base($aln->qseq, $p->qpos, 1);
      } else {
        $temppos = $aln->l_qseq - $p->qpos;
#        ($context[$i], $nt[$i]) = get_context($aln->qseq, $p->qpos, -1);
        ($context[$i], $nt[$i]) = get_2base($aln->qseq, $p->qpos, -1);
      }
      $temppos = -$temppos if $aln->get_tag_values("SECOND_MATE");
      $temphist->{$rg[$i]}->{$temppos}->[$tempq]++;	# position specific
      $temphist->{$rg[$i]}->{$context[$i]}->[$tempq]++;
      if (!$MAXQUAL) {
        $MAXQUAL = $tempq;
      }
      if (!$MINQUAL) {
        $MINQUAL = $tempq;
      }
      if ($tempq > $MAXQUAL) { $MAXQUAL = $tempq; }
      if ($tempq < $MINQUAL) { $MINQUAL = $tempq; }
      $base->{$nt[$i]}++;
    }
    # if we're high coverage, we still need the histograms, but don't use these
    # bases for recalibration - they set this in GATK but unclear it affects us
    return if $MAXCOV && $cov > $MAXCOV;
    return if $MINCOV && $cov < $MINCOV;

    $have_minor = 0;
    foreach $n (qw(G A T C)) {
      $num++ if $base->{$n} > 0;
      if ($base->{$n} > $max) {
        @con = ($n);
        $max = $base->{$n};
      } elsif ($base->{$n} == $max) {
        push @con, $n;
      }
      if ($base->{$n} == 1 ||
          ($base->{$n} > 0 && $cov > 0 && $base->{$n}/$cov < $MINOR_FREQ)
         ) {
        $have_minor = 1;
      }
    } # done with first pass

    # actually only go on if we have minor bases
    if ($have_minor) {
      if (scalar @con) {
        $consensus_base = $con[int(rand(scalar(@con)))];
      } else {
        $consensus_base = ".";
      }

      # randomize the order in case we're not printing all consensus
      $printed_consensus = 0;
      $printed_nonminor = 0;
      @index = ();
      foreach $i (0..$#$pileup) {
        $index[$i] = rand();
      }
      foreach $i (sort { $index[$a] <=> $index[$b] } (0..$#$pileup)) {
        next if !$a[$i];
        $p = $pileup->[$i];
        next if ($p->is_refskip || $p->indel);
        $aln = $a[$i];
        next if !$aln->qual || $aln->qual < $MINMAPQ;
        next if $aln->qscore->[$p->qpos] < $MINQ;
        @data = ();
        push @data, $aln->qscore->[$p->qpos];
        # Consensus or not
        if ($consensus_base eq $nt[$i]) {
          next if ($CONSENSUS_TO_PRINT && $printed_consensus >= $CONSENSUS_TO_PRINT);
          push @data, 1;	# consensus
          $printed_consensus++;
        } else {
          push @data, 0;	# nonconsensus
          if ($cov > 0 && $base->{$nt[$i]} > 1 && $base->{$nt[$i]}/$cov > $MINOR_FREQ) {
            # test coverage for divide by zero (shouldn't happen)
            # we want all singles to get printed regardless of coverage
            # otherwise, limit only when it's not a minor base
            next if ($CONSENSUS_TO_PRINT && $printed_nonminor >= $CONSENSUS_TO_PRINT);
            $printed_nonminor++;
          }
        }
        # "vote" for this nucleotide at this position
        push @data, $base->{$nt[$i]};
        # Coverage at this coordinate
        push @data, $cov;
        # Position - in read
        if ($aln->strand > 0) {
          $temppos = $p->qpos+1;
        } else {
          $temppos = $aln->l_qseq - $p->qpos;
        }
        $temppos = -$temppos if $aln->get_tag_values("SECOND_MATE");
        push @data, $temppos;
        # context - encoded as integer
        push @data, $context_code{$context[$i]};
        if ($DUMP) {
          lock(@DUMPDATA);
          push @DUMPDATA, join ("\t", $seqid, $pos, $nt[$i], $aln->qual, $aln->strand, @data);
        }

        $tempdata->{$rg[$i]} = null if !defined $tempdata->{$rg[$i]};
        $tempdata->{$rg[$i]} = $tempdata->{$rg[$i]}->glue(1,pdl(\@data));
      }
    }
  };	# sub map_consensus

  $sam->clone;
  while ($region = $q->dequeue_nb) {
    last if ($THREADSTATUS[$thread] eq "interrupt");
    print STDERR "$THREADSTATUS[$thread]\n" if $THREADSTATUS[$thread] ne "started";
    $tempdata = {};
    $temphist = {};
    if ($region =~ /:(\d+)-(\d+)/) {
      ($first, $last) = ($1, $2);
    } else {
      ($first, $last) = (0, 0);	# we should never really hit this else clause
    }
if ($thread_debug) { print STDERR "About to do pileup for $region, thread $thread, status $THREADSTATUS[$thread]\n"; }
    $sam->fast_pileup($region, $map_consensus);
if ($thread_debug) { print STDERR "Return from pileup for $region, thread $thread, status $THREADSTATUS[$thread]\n"; }
    foreach $rg (keys %$tempdata) {
      if (!defined $alldata->{$rg}) {
        # actually if there are no read groups we should have defined
        # a default "NULL" read group and set $alldata->{NULL}
        # so the second check here shouldn't really be necessary
        if ($rg eq "NULL" && $RG_LIST[0] ne "NULL") {
          print STDERR "Error, we found some reads not marked by read groups. Ignoring these...\n";
        } else {
          print STDERR "Error, we found some reads marked by a read group ($rg) not in the SAM header. Ignoring these...\n";
        }
        next;
      }
      if ($alldata->{$rg}->at(0,0) == 0) {
        $alldata->{$rg} = $tempdata->{$rg}->copy;
      } else {
        $alldata->{$rg} = $alldata->{$rg}->glue(1, $tempdata->{$rg});
      }
    }
    {
      lock($ALLHIST);
      &add_hist($temphist);
    }
    {
      lock($CALIBRATION_BASES);
      lock($TOTAL_BASES);
      foreach $rg (keys %$tempdata) {
        $CALIBRATION_BASES += $tempdata->{$rg}->getdim(1);
      }
      foreach $rg (keys %$temphist) {
        foreach my $j (0..$#{$temphist->{$rg}->{0}}) {
          $TOTAL_BASES += $temphist->{$rg}->{0}->[$j] if defined $temphist->{$rg}->{0}->[$j];
        }
      }
    }
  }

  # samtools module crashes on thread exit
  # just detach and sleep a long time
  foreach $rg (@RG_LIST) {
    share_pdls($thread . "final" . $rg => $alldata->{$rg}->sever);
  }
  {
    lock @THREADSTATUS;
    $THREADSTATUS[$thread] = "data";
  }
  while (1) {
    sleep 5;
    if ($THREADSTATUS[$thread] eq "goaway") {
      threads->detach;
      last;
    }
  }
  foreach $rg (keys %$alldata) {
    free_pdls($thread . "__" . $rg);
    free_pdls($thread . "final" . $rg);
  }
  while (1) {
    sleep 600;
  }

}

@threadlist = threads->list;
$waiting = 1;
$wait_interval = 60;
if ($verbose) {
  $progress = Term::ProgressBar->new({ name => "Pileups", count => $number_regions, remove => 0, ETA => 'linear' });
  $progress->max_update_rate($wait_interval);
  $progress->minor(0);
  $progress->message("Lacer version $lacer_version running on $bamfile");
  if ($stop_bases > 0) {
    $progress->message("Will stop after finding at least $stop_bases calibration bases. Progress bar won't mean anything...");
  }
  $progress->message("For all read groups: 0/0 (calibration/total; 0.00%) bases processed");
}
if ($DUMP) {
  print STDOUT join ("\t", "# SeqID", "Position", "Base", "MapQ", "Strand", "Quality", "Consensus", "Vote", "Coverage", "Cycle", "ContextCode"), "\n";
}
my $message = "";
my $ml = 0;
while ($waiting) {
  $waiting = 0;
  sleep $wait_interval;
  if ($DUMP && scalar @DUMPDATA) {
    lock(@DUMPDATA);
    print STDOUT join ("\n", @DUMPDATA), "\n";
    @DUMPDATA = ();
  }
  foreach $i (1..$#THREADSTATUS) {
    if ($THREADSTATUS[$i] eq "data") {
      foreach $rg (@RG_LIST) {
        $temppdl = retrieve_pdls($i . "final" . $rg);
        if ($temppdl->at(0,0) > 0) {
          $alldata->{$rg} = $alldata->{$rg}->glue(1, $temppdl);
          $alldata->{$rg}->sever;
        }
        free_pdls($i . "__" . $rg);
      }
      {
        lock @THREADSTATUS;
        $THREADSTATUS[$i] = "goaway";
      }
    }
  }
  foreach $i (@threadlist) {
    $waiting++ if !$i->is_detached;
  }
  if ($verbose) {
    $pending = $queue->pending();
    if (defined($pending)) {
      $progress->update($number_regions - $pending - $num_threads);
    }
    if ($TOTAL_BASES > 0) {
      $message = "\b\rFor all read groups: $CALIBRATION_BASES/" . $TOTAL_BASES . " (calibration/total; " . sprintf("%.2f", $CALIBRATION_BASES/$TOTAL_BASES * 100) . "%) bases processed";
      $ml = length($message) - $ml - 1;
      $progress->message($message);
    }
  }
  if ($thread_debug) {
    print STDERR "Stop Bases: $stop_bases\n";
    print STDERR "Calibration Bases: $CALIBRATION_BASES\n";
    print STDERR "Threads:\n";
    foreach $i (1..$#THREADSTATUS) {
      print STDERR "  Thread $i -- $THREADSTATUS[$i]\n";
    }
    print STDERR "Queue:\n";
    if ($pending) {
      print STDERR "  Pending: $pending; Total: $number_regions\n";
      print STDERR "  Next in queue: ", $queue->peek(), "\n";
      print STDERR "  Last in queue: ", $queue->peek(-1), "\n";
    }
  }
  if ($stop_bases > 0 && $CALIBRATION_BASES > $stop_bases) {
    &quit_pileups;
  }
}

# we should be done now
if (!defined $ALLHIST || ref($ALLHIST) ne "HASH") {
  print STDERR "Uh oh, we don't seem to have any data\n";
  exit;
}
if ($verbose) {
  $progress->update($number_regions);
  print STDERR "\nBeginning SVD-based recalibration...\n";
}

#
# do the recalibration
#
# first clean all the histograms
if ($DUMP) {
  foreach $j ($MINQUAL..$MAXQUAL) {
    foreach $rg (sort keys %$ALLHIST) {
      @out = ("# Quality", "Read Group");
      foreach $covariate (sort keys %{$ALLHIST->{$rg}}) {
        push @out, $covariate;
      }
      print join ("\t", @out), "\n";
      last;
    }
  }
}
foreach $j ($MINQUAL..$MAXQUAL) {
  foreach $rg (sort keys %$ALLHIST) {
    @out = ($j, $rg);
    foreach $covariate (sort keys %{$ALLHIST->{$rg}}) {
      if (!defined $ALLHIST->{$rg}->{$covariate}->[$j] ||
          $ALLHIST->{$rg}->{$covariate}->[$j] < 0) {
        $ALLHIST->{$rg}->{$covariate}->[$j] = 0;
      }
      if ($DUMP) {
        push @out, $ALLHIST->{$rg}->{$covariate}->[$j];
      }
    }
    if ($DUMP) {
      print STDOUT join ("\t", @out), "\n";
    }
  }
}

my $pca1;
foreach $rg (keys %$alldata) {
  # sanity check to make sure we have data
  if (ref($alldata->{$rg}) ne "PDL" ||
      $alldata->{$rg}->isnull ||
      $alldata->{$rg}->at(0,0) == 0) {
    print STDERR "Uh oh, seem to have no data for read group $rg\n";
    next;
  }
  # first the overall recalibration
  ($matrix, $last_count) = make_matrix($alldata->{$rg}, $svd_bin);
  # make the svd_hist here anyway so we can count bases even if we don't have
  # enough data to recalibrate
# $svd_hist->{$rg} = long(@{$ALLHIST->{$rg}->{0}}[$MINQUAL..$MAXQUAL]);
# getting some odd type problem converting to pdl - so make the whole array int?
  $svd_hist->{$rg} = pdl(to_int(@{$ALLHIST->{$rg}->{0}}[$MINQUAL..$MAXQUAL]));
  if ($verbose == 2) {
    $mu->record("After generating matrix");
    print STDOUT "# Initial base data size for read group $rg: ", $alldata->{$rg}->shape, "\n";
  }
  if (!$matrix->isnull) {
    ($recalibrated->{$rg}, $pca1) = quality_svd($matrix, $svd_hist->{$rg}, $svd_bin, $last_count, $min_span);
    if ($verbose) {
      print STDOUT "# SVD fit (overall recalibration): $pca1\n";
    }
    $mu->record("After quality_svd");
  }

  if ($do_covariates) {
    # now covariates - context
    foreach $covariate (list($alldata->{$rg}->(5,)->uniq)) {
#      next if !$covariate;
      $ix = which($alldata->{$rg}->(5,)->flat == $covariate);
      ($matrix, $last_count) = make_matrix($alldata->{$rg}->(:,$ix), $covariate_binsize);
      if (!$matrix->isnull) {
        if (!$gatk) {
          print OUT "# Covariate: $code_context{$covariate} ($covariate)\n";
        }
        $covariate_hist->{$rg}->{$code_context{$covariate}} = pdl(to_int(@{$ALLHIST->{$rg}->{$code_context{$covariate}}}[$MINQUAL..$MAXQUAL]));
        ($covariate_recalibrated->{$rg}->{$code_context{$covariate}}, $pca1) = quality_svd($matrix, $covariate_hist->{$rg}->{$code_context{$covariate}}, $covariate_binsize, $last_count, $min_span);
        if ($verbose) {
          print STDOUT "# SVD fit $code_context{$covariate} ($covariate): $pca1\n";
        }
      }
    }
    $mu->record("After context covariates");

    # now covariates - position
    foreach $covariate (list($alldata->{$rg}->(4,)->uniq)) {
      $ix = which($alldata->{$rg}->(4,)->flat == $covariate);
      ($matrix, $last_count) = make_matrix($alldata->{$rg}->(:,$ix), $covariate_binsize);
      if (!$matrix->isnull) {
        if (!$gatk) {
          print OUT "# Covariate: $covariate\n";
        }
        $covariate_hist->{$rg}->{$covariate} = pdl(to_int(@{$ALLHIST->{$rg}->{$covariate}}[$MINQUAL..$MAXQUAL]));
        ($covariate_recalibrated->{$rg}->{$covariate}, $pca1) = quality_svd($matrix, $covariate_hist->{$rg}->{$covariate}, $covariate_binsize, $last_count, $min_span);
        if ($verbose) {
          print STDOUT "# SVD fit (Position $covariate): $pca1\n";
        }
      }
    }
    $mu->record("After position covariates");
  }
}

#
# output GATK formatted file
#
# first make sure we have something to output
@table = ();
$max_bases = 0;
foreach $rg (keys %$alldata) {
  if (defined $svd_hist->{$rg} &&
      defined $recalibrated->{$rg} &&
      ref($recalibrated->{$rg}) eq "PDL" && 
      ref($svd_hist->{$rg}) eq "PDL") {
    push @table, $rg;
    $max_bases = sum($svd_hist->{$rg}) if sum($svd_hist->{$rg}) > $max_bases;
    if ($verbose == 2) {
      print STDERR "Read group $rg recalibrated on ", $alldata->{$rg}->getdim(1), " nonconsensus and matched consensus bases out of ", sum($svd_hist->{$rg}), " total bases.\n";
    }
  } else {
    print STDERR "== Error ==\n";
    print STDERR "Didn't get recalibration results. This is most likely because we don't have enough data. We need usually at least 5 rows to do an SVD.\n";
    if (defined $svd_hist->{$rg} && ref($svd_hist->{$rg}) eq "PDL") {
      print STDERR " -- For readgroup $rg, we looked at ", sum($svd_hist->{$rg}), " bases total.\n";
      print STDERR "    We identified ", $alldata->{$rg}->getdim(1), " nonconsensus and matched consensus bases to\n";
      print STDERR "    recalibrate on. With an svdbin size of $svd_bin, this gives ", int($alldata->{$rg}->getdim(1)/$svd_bin), " rows for an SVD.\n";
    }
    print STDERR "\n";
    print STDERR "If you don't have enough data, then you can try expanding the region parameter\n";
    print STDERR "or reducing the -svdbin parameter. If you don't see read group specific\n";
    print STDERR "information above, then we probably have a problem with your bam file.\n";
    print STDERR "If you're sure your bam file is ok then please contact the authors.\n";
    print STDERR "Also, if you do have enough data (5 rows) for at least one of the read groups\n";
    print STDERR "above and you are still getting this error message, please contact the authors.\n";
  }
}
exit if !scalar(@table);
$max_bases = length($max_bases);

if ($gatk) {
  &gatk_header($svd_hist, $MAXQUAL, $do_covariates, $max_bases);

  @table = ("ReadGroup    EventType  EmpiricalQuality  EstimatedQReported  Observations  Errors   \n");
  foreach $rg (keys %$recalibrated) {
    next if !defined $svd_hist->{$rg};
    next if ref($recalibrated->{$rg}) ne "PDL";
    next if ref($svd_hist->{$rg}) ne "PDL";
    push @table, gatk_table0($rginfo->{$rg}->{$rgfield}, pdl($MINQUAL..$MAXQUAL), $recalibrated->{$rg}, $svd_hist->{$rg}, $max_bases);
  }
  $left_fields = { 0 => 1,	# ReadGroup
                   1 => 1,	# EventType
                   2 => 0,	# EmpiricalQuality
                   3 => 0,	# EstimatedQReported
                   4 => 0,	# Observations
                   5 => 0 };	# Errors
  &clean_gatk_table($left_fields, \@table);
  # subtract 1 from number of lines in table - we added the header
  print OUT '#:GATKTable:6:' . (scalar(@table)-1) . ':%s:%s:%.4f:%.4f:%d:%.2f:;', "\n";
  print OUT '#:GATKTable:RecalTable0:', "\n";
  print OUT @table;
  print OUT "\n";

  @table = ("ReadGroup    QualityScore  EventType  EmpiricalQuality  Observations  Errors   \n");
  foreach $rg (keys %$recalibrated) {
    next if !defined $svd_hist->{$rg};
    next if ref($recalibrated->{$rg}) ne "PDL";
    next if ref($svd_hist->{$rg}) ne "PDL";
    push @table, gatk_table1($rginfo->{$rg}->{$rgfield}, pdl($MINQUAL..$MAXQUAL), $recalibrated->{$rg}, $svd_hist->{$rg}, $max_bases);
  }
  $left_fields = { 0 => 1,	# ReadGroup
                   1 => 0,	# QualityScore
                   2 => 1,	# EventType
                   3 => 0,	# EmpiricalQuality
                   4 => 0,	# Observations
                   5 => 0 };	# Errors
  &clean_gatk_table($left_fields, \@table);
  # subtract 1 from number of lines in table - we added the header
  print OUT '#:GATKTable:6:' . (scalar(@table)-1) . ':%s:%s:%s:%.4f:%d:%.2f:;', "\n";
  print OUT '#:GATKTable:RecalTable1:', "\n";
  print OUT @table;
  print OUT "\n";

  @table = ("ReadGroup    QualityScore  CovariateValue  CovariateName  EventType  EmpiricalQuality  Observations  Errors  \n");
  foreach $rg (keys %$covariate_recalibrated) {
    if ($do_covariates) {
      foreach $covariate (sort keys %{$covariate_recalibrated->{$rg}}) {
        next if ref($covariate_recalibrated->{$rg}->{$covariate}) ne "PDL";
        next if !defined $covariate_hist->{$rg}->{$covariate};
        next if ref($covariate_hist->{$rg}->{$covariate}) ne "PDL";
        next if $covariate =~ /\./;
        push @table, gatk_table2($rginfo->{$rg}->{$rgfield}, pdl($MINQUAL..$MAXQUAL), $covariate_recalibrated->{$rg}->{$covariate}, $covariate_hist->{$rg}->{$covariate}, $covariate, $max_bases);
      }
    }
  }
  $left_fields = { 0 => 1,	# ReadGroup
                   1 => 0,	# QualityScore
                   2 => 1,	# CovariateValue
                   3 => 1,	# CovariateName
                   4 => 1,	# EventType
                   5 => 0,	# EmpiricalQuality
                   6 => 0,	# Observations
                   7 => 0 };	# Errors
  &clean_gatk_table($left_fields, \@table);
  # subtract 1 from number of lines in table - we added the header
  print OUT '#:GATKTable:8:' . (scalar(@table)-1) . ':%s:%s:%s:%s:%s:%.4f:%d:%.2f:;', "\n";
  print OUT '#:GATKTable:RecalTable2:', "\n";
  print OUT @table;
  print OUT "\n";

}

if ($verbose == 2) {
  $mu->record("Finished");
  @table = split /\n/, $mu->report();
  foreach $i (0..$#table) {
    print STDOUT "# ", $table[$i], "\n";
  }
  print STDOUT "# Total run time ", (time() - $starttime), " seconds.\n";
}

close OUT;

sub get_context {
  my ($sequence, $position, $strand) = @_;
  # position should be 0-based
  # nt should not be reverse complemented - we need to keep track of polymorphism
  # we are giving left and right bases here
  my $left;
  my $right;
  my $temp;
  my $l;
  my $base;

  $l = length($sequence);
  if ($position > 0) {
    ($left, $base) = unpack("x" . ($position - 1) . "A1A1", $sequence);
  } else {
    $left = ".";
    $base = unpack("A1", $sequence);
  }
  if ($position == $l-1) {
    $right = ".";
  } else {
    $right = unpack("x" . ($position + 1) . "A1", $sequence);
  }
  if ($strand < 0) {
    $temp = $left;
    $left = $right;
    $right = $temp;
    $left =~ tr/GATC/CTAG/;
    $right =~ tr/GATC/CTAG/;
  }
  $left =~ s/[^GATC.]/./g;
  $right =~ s/[^GATC.]/./g;
  return ($left.$right, $base);
}

sub get_contextpre {
  my ($sequence, $position, $strand) = @_;
  # position should be 0-based
  # nt should not be reverse complemented - we need to keep track of polymorphism
  # we are giving the two bases before the indicated base here
  my $pre1 = "";
  my $pre2 = "";
  my $l = length($sequence);
  my $base;
  if ($strand > 0) {
    if ($position > 1) {
      ($pre2, $pre1, $base) = unpack("x" . ($position-2) . "A1A1A1", $sequence);
    } elsif ($position == 1) {
      $pre2 = ".";
      ($pre1, $base) = unpack("x" . ($position-1) . "A1A1", $sequence);
    } elsif ($position == 0) {
      $pre2 = ".";
      $pre1 = ".";
      $base = unpack("A1", $sequence);
    }
  } else {
    if ($position < $l-2) {
      if ($position > 0) {
        ($base, $pre1, $pre2) = unpack("x" . ($position) . "A1A1A1", $sequence);
      } else {
        ($base, $pre1, $pre2) = unpack("A1A1A1", $sequence);
      }
    } elsif ($position == $l-2) {
      ($base, $pre1) = unpack("x" . ($position) . "A1A1", $sequence);
      $pre2 = ".";
    } elsif ($position == $l-1) {
      $pre1 = ".";
      $pre2 = ".";
      $base = unpack("x" . ($position) . "A1", $sequence);
    }
    $pre1 =~ tr/GATC/CTAG/;
    $pre2 =~ tr/GATC/CTAG/;
  }
  $pre1 =~ s/[^GATC.]/./g;
  $pre2 =~ s/[^GATC.]/./g;
  return ($pre2.$pre1, $base);
}

sub get_2base {
  my ($sequence, $position, $strand) = @_;
  # position should be 0-based
  # nt should not be reverse complemented - we need to keep track of polymorphism
  # we are giving the preceding and the indicated base here
  my $preceding;
  my $indicated;
  my $base;	# indicated is actually the same as base...unless mapped to
		# negative strand. We want to give positive strand base
		# (which is returned by samtools) but context requires revcomp
  my $temp;
  my $l;

  $l = length($sequence);
  if ($strand > 0) {
    if ($position > 0) {
      ($preceding, $indicated) = unpack("x" . ($position - 1) . "A1A1", $sequence);
    } else {
      $preceding = ".";
      $indicated = unpack("A1", $sequence);
    }
    $base = $indicated;
  } else {
    if ($position >= $l-1) {
      $preceding = ".";
      $indicated = unpack("x" . ($l-1) . "A1", $sequence);
    } elsif ($position == 0) {
      ($indicated, $preceding) = unpack("A1A1", $sequence);
    } else {
      ($indicated, $preceding) = unpack("x" . $position . "A1A1", $sequence);
    }
    $base = $indicated;
    $preceding =~ tr/GATC/CTAG/;
    $indicated =~ tr/GATC/CTAG/;
  }
  $preceding =~ s/[^GATC.]/./g;
  $indicated =~ s/[^GATC.]/./g;
  return ($preceding.$indicated, $base);
}

sub pdlhist {
  my ($q) = @_;
  # give back quality hist respecting minqual and maxqual
  # assume qualities are in the first column
  my $min = minimum($q);
  my $max = maximum($q) + 1;
  # we get bins starting at $min+0.5 and going to $max-0.5
  # quality score of q goes into bin q+0.5 - so need max of quality + 1
  my ($bin, $hist) = hist($q, $min, $max, 1);
  # below are actual qualities that go into bins
  my $minhist = int($bin->at(0));
  my $maxhist = int($bin->at(-1));
  if ($MINQUAL < $minhist) {
    $hist = append(zeroes($minhist - $MINQUAL), $hist);
  } elsif ($MINQUAL > $minhist) {
    $hist = $hist->(($MINQUAL - $minhist):);
  }
  if ($MAXQUAL > $maxhist) {
    $hist = append($hist, zeroes($MAXQUAL - $maxhist));
  } elsif ($MAXQUAL < $maxhist) {
    $hist = $hist->(:($MAXQUAL - $maxhist));
  }
  return $hist/sum($hist);
}

sub quality_svd {
  my ($matrix, $hist, $bin_size, $last_count, $min_span) = @_;

  # some variables
  my $overall_q = $hist / sum($hist);
  my $total_bases = sum($hist);
  my $tosvd;
  my $center;
  my ($u, $s, $v);
  my $fit;
  my $candidates;
  my $tolerance;
  my $sdev;
  my ($c1, $c2);
  my $candix;
  my $correct;
  my $error;
  my $correct_x;
  my $error_x;
  my $last_row;
  my $total_error;
  my $recalibrated;
  my $errorperc;
  my $u0;
  my ($n, $m, $i, $append);

  # center the columns of the matrix, do the svd
  $tosvd = $matrix->(,0:-2)->copy;
  $center = average($tosvd->xchg(0,1));
  $tosvd -= $center;
  $n = $tosvd->dim(0);
  $m = $tosvd->dim(1);
  # somewhere between 2.007 and 2.018 there is a check for $m > $n in the svd. Just copy the rows until we have enough.
  if ($m < $n) {
    $append = $tosvd->copy;
    foreach $i (1..int($n/$m)) {
      $tosvd = $tosvd->glue(1, $append);
    }
  }
  ($u, $s, $v) = svd($tosvd);
  $u = $u->(,0:($m-1));
  $fit = $s->at(0)*$s->at(0)/sum($s*$s);

  # try to find the best min and max dim 0 coordinate
  $candidates = $u->(0,)->(0);
  $tolerance = 0;
  $sdev = sqrt(inner($u->(0,)->flat,$u->(0,)->flat)/$u->getdim(1));
  while ($tolerance <= 2 &&
         ($candidates->getdim(0) < 2 ||
          max($candidates) - min($candidates) < $min_span * (max($u->(0,)) - min($u->(0,)))
        )) {
    $tolerance += 0.1;
    $c1 = which(abs($u->(1,)) < $tolerance * $sdev);
    $c2 = which(abs($u->(2,)) < $tolerance * $sdev);
    if ($c1->getdim(0) > 0 && $c2->getdim(0) > 0) {
      $candix = intersect($c1, $c2);
    } else {
      $candix = pdl(0);
    }
    if ($candix->getdim(0) == 0) {
      $candix = pdl(0);
    }
    $candidates = $u->(0,)->flat->index($candix);
  }
  if ($candidates->getdim(0) < 2) {
    $candidates = $u->(0,)->flat;
  }

  # look at loading direction of highest quality base
  if ($v->at(0,-1) > 0) {
    $correct = $center + max($candidates) * $v->(0,)->flat * $s->at(0);
    $error = $center + min($candidates) * $v->(0,)->flat * $s->at(0);
    $correct_x = max($candidates);
    $error_x = min($candidates);
  } else {
    $correct = $center + min($candidates) * $v->(0,)->flat * $s->at(0);
    $error = $center + max($candidates) * $v->(0,)->flat * $s->at(0);
    $correct_x = min($candidates);
    $error_x = max($candidates);
  }

  # clip to min and max of candidates for error calculation
  $u0 = $u->(0,)->flat->copy;
  $errorperc = $u0->inplace->clip(min($candidates), max($candidates));
  $errorperc = ($errorperc - $correct_x) / ($error_x - $correct_x);
  # calculate error bases, take care of last row
  $total_error = sum($errorperc) * $bin_size;
  $last_row = ($matrix->(,-1)->flat - $center) x $v;
  $last_row /= $s->(0);
  if ($last_row->at(0,0) > max($candidates)) {
    $last_row->(0,0) .= max($candidates);
  } elsif ($last_row->at(0,0) < min($candidates)) {
    $last_row->(0,0) .= min($candidates);
  }
  $total_error += $last_count * ($last_row->at(0,0) - $correct_x) / ($error_x - $correct_x);

  # "yates" correction like GATK
  $error *= $total_error;
  $error->inplace->clip(1,undef);
  $error /= sum($error);
  $correct *= ($total_bases - $total_error);
  $correct->inplace->clip(1,undef);
  $correct /= sum($correct);
  $overall_q *= $total_bases;
  $overall_q->inplace->clip(2, undef);
  $overall_q /= sum($overall_q);

  # the recalibration
  $recalibrated = -10 * log10($total_error/$total_bases * $error / $overall_q);
  if (!$gatk) {
    print OUT "# Initial matrix size ", $matrix->shape, "\n";
    print OUT "# Bin size: $bin_size\n";
    print OUT "# SVD fit: ", $s->at(0) * $s->at(0) / inner($s, $s), "\n";
    print OUT "# Tolerance: $tolerance\n";
    print OUT "# Stdev: $sdev\n";
    print OUT "# Candidates: $candidates\n";
    print OUT "# U1 range: ", (max($u->(0,)) - min($u->(0,))), "\n";
    print OUT join ("\t", "# Original", "Recalibrated", "Histogram", "Correct", "Error"), "\n";
    wcols(pdl($MINQUAL..$MAXQUAL), $recalibrated, $hist, $correct, $error, { COLSEP => "\t" });
  }
  
  return($recalibrated, $s->at(0) * $s->at(0) / inner($s, $s));
}

sub gatk_header {
  my ($hist, $MAXQUAL, $do_covariates, $max_length) = @_;
  # just dummy stuff
  if ($do_covariates) {
    print OUT '#:GATKReport.v1.1:5
#:GATKTable:2:17:%s:%s:;
#:GATKTable:Arguments:Recalibration argument collection values used in this run
Argument                    Value                                                                   
binary_tag_name             null                                                                    
covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate
default_platform            null                                                                    
deletions_default_quality   45                                                                      
force_platform              null                                                                    
indels_context_size         3                                                                       
insertions_default_quality  45                                                                      
low_quality_tail            2                                                                       
maximum_cycle_value         500                                                                     
mismatches_context_size     2                                                                       
mismatches_default_quality  -1                                                                      
no_standard_covs            true                                                                    
quantizing_levels           16                                                                      
recalibration_report        null                                                                    
run_without_dbsnp           false                                                                   
solid_nocall_strategy       THROW_EXCEPTION                                                         
solid_recal_mode            SET_Q_ZERO                                                              

#:GATKTable:3:94:%s:%s:%s:;
#:GATKTable:Quantized:Quality quantization map
';
  } else {
    print OUT '#:GATKReport.v1.1:5
#:GATKTable:2:17:%s:%s:;
#:GATKTable:Arguments:Recalibration argument collection values used in this run
Argument                    Value                                   
binary_tag_name             null                                    
covariate                   ReadGroupCovariate,QualityScoreCovariate
default_platform            null                                    
deletions_default_quality   45                                      
force_platform              null                                    
indels_context_size         3                                       
insertions_default_quality  45                                      
low_quality_tail            2                                       
maximum_cycle_value         500                                     
mismatches_context_size     2                                       
mismatches_default_quality  -1                                      
no_standard_covs            true                                    
quantizing_levels           16                                      
recalibration_report        null                                    
run_without_dbsnp           false                                   
solid_nocall_strategy       THROW_EXCEPTION                         
solid_recal_mode            SET_Q_ZERO                              

#:GATKTable:3:94:%s:%s:%s:;
#:GATKTable:Quantized:Quality quantization map
';
  }

  my $sum = 0;
  my $rg;
  foreach $rg (keys %$hist) {
    $sum += sum($hist->{$rg});
  }
  my $width = 11;
  $width = $max_bases + 2 if $max_bases + 2 > $width;
  print OUT "QualityScore  Count" . (" " x ($width-7)) . "  QuantizedScore\n";
  foreach my $i (0..93) {
    if ($i < $MAXQUAL) {
      printf OUT ("% 12d% *d% 16d\n", $i, $width, 0, 8);
    } elsif ($i > $MAXQUAL) {
      printf OUT ("% 12d% *d% 16d\n", $i, $width, 0, $MAXQUAL);
    } else {
      printf OUT ("% 12d% *d% 16d\n", $i, $width, $sum, $MAXQUAL);
    }
  }
  print OUT "\n";
}

sub gatk_table0 {
  # overall quality
  my ($rg, $qbasis, $empirical, $hist, $max_length) = @_;
  my $error = $hist * (10 ** (-$empirical/10));
  my $reportederror = $hist * (10 ** (-$qbasis/10));
  my @return = ();
  my $width = 11;
  $width = $max_bases + 2 if $max_bases + 2 > $width;

  push @return, sprintf("%-12s M% 21.0f.0000% 20.4f% 14d% *.2f\n", $rg, -10*log(sum($error)/sum($hist)) / log(10), -10*log(sum($reportederror)/sum($hist)) / log(10), sum($hist), $width, sum($error));

  return @return;
}

sub gatk_table1 {
  # break out by quality scores
  my ($rg, $qbasis, $empirical, $hist, $max_length) = @_;
  my $error = $hist * (10 ** (-$empirical/10));
  my @error = list($error);
  my @qbasis = list($qbasis);
  my @empirical = list($empirical);
  my @hist = list($hist);
  my @return = ();
  my $width = 11;
  $width = $max_bases + 1 if $max_bases + 2 > $width;

  foreach $i (0..$#qbasis) {
    next if !$hist[$i];
      push @return, sprintf("%-12s%13d  M% 26.4f% 14d% *.2f\n", $rg, $qbasis[$i], $empirical[$i], $hist[$i], $width, $error[$i]);
  }

  return @return;
}

sub gatk_table2 {
  # break out by covariates
  my ($rg, $qbasis, $empirical, $hist, $covariate, $max_length) = @_;
  my $error = $hist * (10 ** (-$empirical/10));
  my @error = list($error);
  my @qbasis = list($qbasis);
  my @empirical = list($empirical);
  my @hist = list($hist);
  my @return = ();
  my $width = 10;
  $width = $max_bases + 1 if $max_bases + 2 > $width;

  foreach $i (0..$#qbasis) {
    next if !$hist[$i];
    if ($covariate =~ /[GATC]/) {
      push @return, sprintf("%-12s%13d  %-14s  Context        M% 26.4f% 14d% *.2f\n", $rg, $qbasis[$i], $covariate, $empirical[$i], $hist[$i], $width, $error[$i]);
    } else {
      push @return, sprintf("%-12s%13d  %-14s  Cycle          M% 26.4f% 14d% *.2f\n", $rg, $qbasis[$i], $covariate, $empirical[$i], $hist[$i], $width, $error[$i]);
    }
  }

  return @return;
}

sub make_matrix {
  my ($data, $svd_bin) = @_;
  # we're expecting these fields specified by first dimension coordinate
  # 0 - quality
  # 1 - consensus
  # 2 - vote
  # 3 - coverage
  # 4 - position
  # 5 - encoded context
  my $matrix = null;
  my $scale;
  my $score;
  my $ix;
  my $qonly;
  my $j;
  my $jend;
  my $last_count;
  my $max_cov;
  $max_cov = max($data->(3,));
  $scale = length($max_cov);
  $score = ($max_cov-$data->(2,)) * 10**($scale+2) +	# vote (low then high)
           (1-$data->(1,)) * 10**($scale+1) +	# consensus (0 then 1)
           $data->(3,);				# coverage (high then low)
  $ix = qsorti($score->flat);			# sort ascending
  # but now reverse indices - so descending (i.e. order indicated above)
  $qonly = $data->(0,:)->flat->index($ix->(-1:0))->copy;
  for ($j = 0; $j < $data->getdim(1); $j += $svd_bin) {
    if ($j + $svd_bin < $data->getdim(1)) {
      $jend = $j + $svd_bin - 1;
    } else {
      $jend = $data->getdim(1)-1;
    }
    $matrix = $matrix->glue(1, pdlhist($qonly->($j:$jend)));
    $last_count = $jend - $j + 1;
  }
  if ($matrix->getdim(1) < 5) {
    return (null, 0);
  } else {
    return ($matrix->copy, $last_count);
  }
}

sub add_hist {
  my ($new) = @_;	# references to histograms
  my $rg;
  my $cov;
  my $q;
  foreach $rg (keys %$new) {
    if (!defined($ALLHIST->{$rg})) {
      my %anon :shared;
      $ALLHIST->{$rg} = \%anon;
    }
    foreach $cov (keys %{$new->{$rg}}) {
      if (!defined($ALLHIST->{$rg}->{$cov})) {
        my @anon2 :shared;
        $ALLHIST->{$rg}->{$cov} = \@anon2;
      }
      foreach $q (0..$#{$new->{$rg}->{$cov}}) {
        if (defined $new->{$rg}->{$cov}->[$q]) {
          if (!defined $ALLHIST->{$rg}->{$cov}->[$q]) {
            $ALLHIST->{$rg}->{$cov}->[$q] = $new->{$rg}->{$cov}->[$q];
          } else {
            $ALLHIST->{$rg}->{$cov}->[$q] += $new->{$rg}->{$cov}->[$q];
          }
        }
      }
    }
  }
}

sub to_int {
  my (@a) = @_;
  foreach my $i (0..$#a) {
    $a[$i] += 0;
  }
  return @a;
}

sub sigint_handler {
  if (!$ABORT) {
    print STDERR "\nCaught Control-C.  Finishing current pileup windows and recalibrating.\nPress Control-C again to abort immediately.\n\n";
    &quit_pileups;
    $ABORT = 1;
  } else {
    die;
  }
}

sub quit_pileups {
  my $i;
  foreach $i (1..$#THREADSTATUS) {
    $THREADSTATUS[$i] = "interrupt" if $THREADSTATUS[$i] eq "started";
  }
}

sub clean_gatk_table {
  my ($left, $out) = @_;
  # we already have newlines on these
  # mainly deal with spacing issues
  # need 2 spaces between all fields
  # numbers are mostly right-justified
  # text is left-justified
  # headings are left-justified
  # we need to know which fields are left beforehand, changes with table
  my @width;
  my $i;
  my $j;
  my @f;
  my $new;
  foreach $i (0..$#$out) {
    next if $out->[$i] =~ /^#/;
    next if $out->[$i] =~ /^$/;
    chomp $out->[$i];
    @f = split /\s+/, $out->[$i];
    foreach $j (0..$#f) {
      $width[$j] = 0 if !defined $width[$j];
      $width[$j] = length($f[$j]) if $width[$j] < length($f[$j]);
    }
  }

  # first line should be header, always left justified
  @f = split /\s+/, $out->[0];
  $new = sprintf("%-*s", $width[0], $f[0]);
  foreach $j (1..$#f) {
    $new .= "  ";
    $new .= sprintf("%-*s", $width[$j], $f[$j]);
  }
  $new .= "\n";
  $out->[0] = $new;

  foreach $i (1..$#$out) {
    @f = split /\s+/, $out->[$i];
    # field 0 is always text even if it's numbers - read group
    $new = sprintf("%-*s", $width[0], $f[0]);
    foreach $j (1..$#f) {
      $new .= "  ";	# spacer
      if ($left->{$j}) {
        $new .= sprintf("%-*s", $width[$j], $f[$j]);
      } else {
        $new .= sprintf("%*s", $width[$j], $f[$j]);
      }
    }
    $new .= "\n";
    $out->[$i] = $new;
    $new = "";
  }
}


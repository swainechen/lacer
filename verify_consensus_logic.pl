#!/usr/bin/perl
use strict;
use warnings;
use Test::More;

# Mocking the environment for map_consensus worker logic
my $MINMAPQ = 30;
my $MINQ = 6;
my $USE_READGROUPS = 1;
my %VALID_RG = ( "RG1" => 1 );
my $MINOR_FREQ = 0.05;
my $CONSENSUS_TO_PRINT = 1;

# Mock Alignment object
package MockAln;
sub new {
    my ($class, %args) = @_;
    return bless \%args, $class;
}
sub qual { $_[0]->{qual} }
sub _qscore { $_[0]->{qscore} }
sub aux { $_[0]->{aux} }
sub strand { 1 }
sub l_qseq { length($_[0]->{qseq}) }
sub qseq { $_[0]->{qseq} }
sub get_tag_values { 0 }
sub qscore {
    my ($self) = @_;
    my @scores = map { ord($_) } split //, $self->{qscore};
    return \@scores;
}

package MockPileup;
sub new {
    my ($class, %args) = @_;
    return bless \%args, $class;
}
sub is_refskip { 0 }
sub indel { 0 }
sub alignment { $_[0]->{aln} }
sub qpos { $_[0]->{qpos} }

package main;

sub test_logic {
    my ($pileup_records) = @_;

    my @a = (0) x scalar(@$pileup_records);
    my @rg;
    my @pos;
    my @context;
    my @nt;
    my $cov = 0;
    my $base = { "G" => 0, "A" => 0, "T" => 0, "C" => 0 };
    my $temphist = {};
    my $MAXQUAL = 0;
    my $MINQUAL = 0;

    # First Pass (Simulating updated lacer.pl logic)
    for (my $i = scalar(@$pileup_records)-1; $i >= 0; $i--) {
        my $p = $pileup_records->[$i];

        my $this_aln = $p->alignment;
        next if !$this_aln->qual || $this_aln->qual < $MINMAPQ;

        my $q_scores = $this_aln->_qscore;
        my $tempq = ord(substr($q_scores, $p->qpos, 1));
        next if $tempq < $MINQ || $tempq > 93;

        my $this_rg = "NULL";
        if ($USE_READGROUPS) {
            $this_rg = $this_aln->aux;
            if ($this_rg =~ /RG:Z:(\S+)/) {
                $this_rg = $1;
            }
        }

        # Security validation
        if (!exists $VALID_RG{$this_rg}) {
            # warn "Skipping invalid RG $this_rg\n";
            next;
        }

        # Simplified coordinate/base logic
        my $this_nt = substr($this_aln->qseq, $p->qpos, 1);
        my $this_context = "..";
        my $this_temppos = $p->qpos + 1;

        # State updates (deferred)
        $a[$i] = $this_aln;
        $rg[$i] = $this_rg;
        $pos[$i] = $this_temppos;
        $context[$i] = $this_context;
        $nt[$i] = $this_nt;
        $cov++;

        $base->{$this_nt}++;
    }

    # Second Pass (Simulating lacer.pl logic)
    my $processed_count = 0;
    my $skipped_count = 0;

    # Simple consensus for testing
    my $consensus_base = "A";

    foreach my $i (0..$#$pileup_records) {
        if (!$a[$i]) {
            $skipped_count++;
            next;
        }
        $processed_count++;
    }

    return ($processed_count, $skipped_count);
}

# --- Test Cases ---

# 1. Valid record
my $aln1 = MockAln->new(qual => 40, qscore => chr(30), aux => "RG:Z:RG1", qseq => "A");
my $p1 = MockPileup->new(aln => $aln1, qpos => 0);

my ($proc, $skip) = test_logic([$p1]);
is($proc, 1, "Processed valid record");
is($skip, 0, "No records skipped");

# 2. Invalid Read Group record (Security Bypass Test)
my $aln2 = MockAln->new(qual => 40, qscore => chr(30), aux => "RG:Z:EVIL", qseq => "G");
my $p2 = MockPileup->new(aln => $aln2, qpos => 0);

($proc, $skip) = test_logic([$p2]);
is($proc, 0, "Did NOT process invalid RG record");
is($skip, 1, "Skipped invalid RG record");

# 3. Mixed records
($proc, $skip) = test_logic([$p1, $p2]);
is($proc, 1, "Processed only valid record in mixed set");
is($skip, 1, "Skipped invalid record in mixed set");

done_testing();

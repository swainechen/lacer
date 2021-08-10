# lacer
Lacer: Accurate Base Quality Score Recalibration using Linear Algebra  
Version: 0.426

Lacepr: A fast replacement to rewrite recalibrated base quality scores in bam files  
Version: 0.2

Lacer takes a BAM file and produces a recalibration file similar to GATK's BaseRecalibrator. The algorithm is distinct from GATK's, however, and is usable on non-human data without requirement to already have a set of known variants. The resulting recalibration file can be used in GATK's PrintReads to perform a recalibration.

Lacepr takes a BAM file and recalibration file and produces a new BAM file with recalibrated quality scores. The algorithm is the same as that used in GATK. It's much (~2-3x) faster but does less error checking. It's also a bit forgiving with ReadGroup issues.

Installation / Requirements
---------------------------

Lacer is a perl script.  It has only been tested on Perl version 5.10 and later.

It depends on the following Perl modules available in a standard Perl installation or on CPAN:
- Bio::DB::Sam
- Data::Dumper
- Getopt::Long
- Memory::Usage
- PDL
- PDL::MatrixOps
- PDL::NiceSlice
- PDL::Parallel::threads
- Term::ProgressBar
- Thread::Queue
- threads
- threads::shared
-
To check the availability of these modules, you can run:
```
perl -e 'use Bio::DB::Sam'
```

replacing "Bio::DB::Sam" with any of the modules from the above list.  If there is no error, that module is available.

For Lacepr, one needs the samtools and htslib libraries. It has only been tested with samtools-1.1 and htslib-1.1 and later. To compile, first set the location of the samtools and htslib libraries, then execute the following command in the lacepr subdirectory:
```
SAMTOOLS=path-to-samtools
HTSLIB=path-to-htslib
gcc -I$SAMTOOLS -I$HTSLIB lacepr.c -L$SAMTOOLS -L$HTSLIB -lbam -l:libhts.a -lz -lpthread -lm -o lacepr
```

Usage
-----
A help screen describing options to Lacer is displayed if you run the script without any arguments:
```
./lacer.pl
```

An indexed bam file is needed.  If your bam file is named "input.bam", please run:
```
samtools index input.bam
```

or the equivalent prior to running Lacer.

The reference fasta file to which your reads were aligned is also required.  If this is named "reference.fasta", then the basic usage is:
```
./lacer.pl -bam input.bam -reference reference.fasta -outfile recal.txt
```

The output file, recal.txt, can then be used in the GATK workflow:
```
java -Xmx2g -jar GenomeAnalysisTK.jar \
     -R reference.fasta \
     -T PrintReads \
     -o recal.bam \
     -I input.bam \
     -BQSR recal.txt
```

Alternatively, the provided lacepr command can be used:
```
./lacepr --bam input.bam --recal recal.txt --out recal.bam
```

Note: The format of the recalibration table required by GATK has changed frequently in recent versions.  Depending on the version of GATK you are using, you may need to specify additional parameters to Lacer to produce the appropriate format file.
For GATK version 2.7:
```
./lacer.pl -bam input.bam -reference reference.fasta -rgfield ID -outfile recal.txt
```

For GATK version 2.8 or later:
```
./lacer.pl -bam input.bam -reference reference.fasta -rgfield PU -outfile recal.txt
```

If you are using `lacepr` to print out recalibrated reads, simply use the defaults for Lacer. There are also options for forcing the use of read groups if needed.

Known Issues
------------
When running Lacer, the following warning will appear:
Subroutine PDL::CLONE_SKIP redefined at /mnt/software/lib/perl5/5.10.1/PDL/Parallel/threads.pm line 39.

This warning can be safely ignored.

Contact/citation/information
-------------------
slchen@gis.a-star.edu.sg, swainechen@gmail.com
Preprint: https://doi.org/10.1101/130732
Citation: Chung JC and Chen SL. 2017. "Lacer: Accurate Base Quality Score Recalibration For Improving Variant Calling From Next-Generation Sequencing Data In Any Organism." bioRxiv doi: https://doi.org/10.1101/130732

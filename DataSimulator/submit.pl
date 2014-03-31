#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
Getopt::Long::config("no_ignore_case");
use POSIX qw/ceil/;
use Scalar::Util qw(looks_like_number);


#set dumper for dumping nested data structures, usage: $dump->dumpValue(\@array);
require Dumpvalue;
my $dump = new Dumpvalue::
			    tick        => q("),
			    compactDump => 1,  # comment these two lines out
			    veryCompact => 1,  # if you want a bigger dump
			;
			
my $USAGE = "
USAGE: perl $0 [OPTIONS] <RNArobo binary> <settings file> <in folder> <out folder>

  A script to run RNArobo/RNAbob.

Available options:
  -f 'bflags' options for the binary
          
  -qsub       queue the tests to qsub
  
  -h, --help  show this help

Example usage:
  perl $0 -f '-c' ../rnarobo in/settings.txt in/ out/
  perl $0 -qsub -f '-c' ../rnarobo inSVK/settingsWo.txt inSVK/ outSVK/
  perl $0 -qsub -f '-descr' bin/rnamotif  inSVK/settingsRNAMotif.txt inSVK/ outSVK/
";

my $bflags = '';
my $toQsub = 0;
my $showHelp = 0;
GetOptions( 'f=s' => \$bflags,
            'qsub' => \$toQsub,
            'help' => \$showHelp ) 
      || die $USAGE;
      
if ($showHelp) { die $USAGE; }

#parse command line arguments
#print join ("|", @ARGV), "\n";
my $binary;
my $settingsFileName;
my $testDir;
my $outDir;
if (scalar @ARGV == 4){
    $binary = shift @ARGV;
    $settingsFileName = shift @ARGV;
    $testDir = shift @ARGV;
    $outDir = shift @ARGV;
    
    $testDir .= '/' if $testDir !~ /\/$/;
    $outDir .= '/' if $outDir !~ /\/$/;
} else {
    die $USAGE;
}

################################################################################

#print "$settingsFileName\n";
open FILE, "<$settingsFileName" or die $!;

my $line;
while($line = <FILE>){
  $line =~ s/^\s*//;
  $line =~ s/\s*$//;
  next if $line =~ /^\#/;
  
  my ($des, $db) = split ' ', $line;
  #print "$des $db\n";
  
  my $out = "";
  $out .= $1 if $des =~ /\/?(.+)\.[^\.]+$/;
  #print "$out\n";
  
  print "Running $binary $bflags on $testDir$des $testDir$db > $outDir$out\n";
  if($toQsub){
    `qsub -q all.q -R y -pe parallel 6 -l h_vmem=10G -l h_rt=05:00:00 -S /bin/bash some_shell_script.sh param1 param2`;
  } else {
    `/usr/bin/time -p $binary $bflags $testDir$des $testDir$db > $outDir$out 2>> $outDir$out &`;
  }
}

close FILE;



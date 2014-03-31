use strict;
use Getopt::Long;

my $USAGE = "
USAGE: perl $0 [OPTIONS] [-i <input_files_prefix>] [-b <begin num>] [-e <end_num>]\n
";

#parse command line arguments
#die $USAGE unless @ARGV>=1;

#parse commandline parameters - input filename, "insert ampersand", ...
my $filenamePrefix = "";
my $beginFNum = 0;
my $endFNum = 0;
my $showHelp = 0;
GetOptions( 'i=s' => \$filenamePrefix,
            'b:i' => \$beginFNum,
            'e:i' => \$endFNum,
            'help' => \$showHelp )
      || die $USAGE;
      
if ($showHelp) { die $USAGE; }

#print "$beginFNum, $endFNum\n";

my @VR_type_ok;
my @VP_type_ok;
my @PR_type_ok;
my @PP_type_ok;
my @VR_type_wr;
my @VP_type_wr;
my @PR_type_wr;
my @PP_type_wr;
my %VR_length_ok;
my %VP_length_ok;
my %PR_length_ok;
my %PP_length_ok;
my %VR_length_wr;
my %VP_length_wr;
my %PR_length_wr;
my %PP_length_wr;

for (my $fileNum = $beginFNum; $fileNum <= $endFNum; $fileNum++) {
 	my $filename = "$filenamePrefix$fileNum.txt";
 	#print $filename, "\n";
 	my $file;
    unless (open $file, "<$filename") {
        warn "Cannot open $filename: $!" unless !$filename;
        next;
    }

    my $active = 0;
    my $counter = 1;
    while(my $line = <$file> ) {
        chomp $line;
        if ($line) {
            #skip this line if it doesnt look like new occurrence
            next unless $line =~ /^([V,P])([R,P])\s(\d+)\s:\s+(\d+)\s+(\d+)/;
            my $inferenceType = $1;
            my $measure = $2;
            my $unit = $3;
            my $correct = $4;
            my $wrong = $5;
            #print "$inferenceType $measure $unit $correct $wrong\n";
            my $pref = $inferenceType.$measure;
            if ($unit < 10){
                if ($pref eq "VR"){
                    $VR_type_ok[$unit] += $correct;
                    $VR_type_wr[$unit] += $wrong;
                } elsif ($pref eq "VP") {
                    $VP_type_ok[$unit] += $correct;
                    $VP_type_wr[$unit] += $wrong;
                } elsif ($pref eq "PR") {
                    $PR_type_ok[$unit] += $correct;
                    $PR_type_wr[$unit] += $wrong;
                } elsif ($pref eq "PP") {
                    $PP_type_ok[$unit] += $correct;
                    $PP_type_wr[$unit] += $wrong;
                }
            } else {
                if ($pref eq "VR"){
                    $VR_length_ok{$unit} += $correct;
                    $VR_length_wr{$unit} += $wrong;
                } elsif ($pref eq "VP") {
                    $VP_length_ok{$unit} += $correct;
                    $VP_length_wr{$unit} += $wrong;
                } elsif ($pref eq "PR") {
                    $PR_length_ok{$unit} += $correct;
                    $PR_length_wr{$unit} += $wrong;
                } elsif ($pref eq "PP") {
                    $PP_length_ok{$unit} += $correct;
                    $PP_length_wr{$unit} += $wrong;
                }
            }
        }
    }
    
    close $file;
}

print "Viterbi:\n";
print "recall by type:\n";
for (my $i=0; $i < 7; $i++){
    my $o = $VR_type_ok[$i];
    my $w = $VR_type_wr[$i];
    print "$i $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}
print "precision by type:\n";
for (my $i=0; $i < 7; $i++){
    my $o = $VP_type_ok[$i];
    my $w = $VP_type_wr[$i];
    print "$i $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}
print "recall by length:\n";
for (sort {$a<=>$b} keys %VR_length_ok){
    my $o = $VR_length_ok{$_};
    my $w = $VR_length_wr{$_};
    print "$_ $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}
print "precision by length:\n";
for (sort {$a<=>$b} keys %VP_length_ok){
    my $o = $VP_length_ok{$_};
    my $w = $VP_length_wr{$_};
    print "$_ $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}

print "\n\nPosterior:\n";
print "recall by type:\n";
for (my $i=0; $i < 7; $i++){
    my $o = $PR_type_ok[$i];
    my $w = $PR_type_wr[$i];
    print "$i $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}
print "precision by type:\n";
for (my $i=0; $i < 7; $i++){
    my $o = $PP_type_ok[$i];
    my $w = $PP_type_wr[$i];
    print "$i $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}
print "recall by length:\n";
for (sort {$a<=>$b} keys %PR_length_ok){
    my $o = $PR_length_ok{$_};
    my $w = $PR_length_wr{$_};
    print "$_ $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}
print "precision by length:\n";
for (sort {$a<=>$b} keys %PP_length_ok){
    my $o = $PP_length_ok{$_};
    my $w = $PP_length_wr{$_};
    print "$_ $o $w ", $o+$w==0 ? 0 : $o/($o+$w)*100, "\n";
}




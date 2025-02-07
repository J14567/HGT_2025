#The output of this file is a table with two column. One column being the position of the read in the genome and the other column the position of specific reads (when existing) on the alternative genome. 
#The obtained file can be analyze on Analysis_original_file.R or alternative of this file. 

use strict;
use warnings;


my $infile=$ARGV[0]; #nb d'arguments
my $sufix1=$ARGV[1];
my $sufix2=$ARGV[2];
$infile=~/(.+).sam/; #récupérer le nom du fichier 

my $prefix=$1;
#my $prefix=$ARGV[];
print "$infile\t$prefix\n";
my $outfile1=$prefix."_$sufix1.sam";  #$ means the value is a scalar
my $outfile2=$prefix."_$sufix2.sam";
my $outfile3=$prefix."_$sufix1"."_$sufix2.sam";
my $outfile4=$prefix."_unmapped.sam";
my $outfile5=$prefix."_Plasmid.sam";



my $line;
my $CMI;
my $O1;
my $O2;
my $O3;
my $O4;
my $O5;

my @b;
open($CMI,'<',$infile) or die("open: $!");  
# open: open(filehandle, mode, filename) avec > write, < read
# die: file-opening failure
open($O1,'>',$outfile1) or die("open: $!");
open($O2,'>',$outfile2) or die("open: $!");
open($O3,'>',$outfile3) or die("open: $!");
open($O4,'>',$outfile4) or die("open: $!");
open($O5,'>',$outfile5) or die("open: $!");

my $target=0;
my $score;
my $score_alternative;
my $secondmatch;
my $nb1=0;
my $nb2=0;
my $nb3=0;
my $nbp=0;
my $unmapped=0;
my $counter=0;
my  $secondpos;
my $multiplehits;

while(defined ($line = <$CMI>))
	{
        
    if ($counter % 1000000==0) 
    {print "$counter\n";}
    $counter++;
	chomp $line;
    if (substr($line,0,1)eq '@')
        {
        next;
        }
    $secondpos=-1;
	@b=split(/\t/,$line);
	#print "$line\n";
        
        #read score
        $line=~/AS:i:(\d+)/;
        $score=$1;
        
        #read score of alternative mapping
        $line=~/XS:i:(\d+)/;
        $score_alternative=$1;
        
        if ($b[2] eq "*" || $score==0)# eq: check equality between two strings 
        {
            $unmapped=$unmapped+1;
            $target=0;
            print $O4 "$line\n";
            next;
        }
        if ($b[2] eq "$sufix1")
            {
            $target=1;
            }
        if ($b[2] eq "$sufix2")
            {
            $target=2;
            }
        if ($b[2] eq "Plasmid")
            {
            $target=4;
            }
        
        #if there is a lesser mapping but relavant enough to be noted we have
        if ($line=~/XA:Z:(.+),.(\d+),.*,\d+;/ && $score_alternative>$score-30)
        {
        $secondmatch=$1;
            if ($secondmatch ne $b[2]) #if it matches on the other genome store the ratio
            #on peut trouver dans les deux génomes 
        #4 options: mapper gén1, mapper gén2, 
            {
             $secondpos=$2; 
            }
            
        #if there are more than two alternative hits
        $multiplehits=0;
            if ($line=~/XA:Z:.+,.\d+,.*,\d+;.+,.\d+,.*,\d+;/)
                {
                    $multiplehits=1;
                    $target=3; # match plusieurs fois dans le génome = poubelle 
                }
            
        
        }
        
        
        #if the alternative mapping is as good
        if ($score_alternative==$score && $score>0) 
        {
            $target=3;
        }
        
        
       
        
        if ($target==1)
        {
            print $O1 "$b[3]\t$secondpos\n";
            $nb1++;
        }
        if ($target==2)
        {
            print $O2 "$b[3]\t$secondpos\n"; 
            $nb2++;
        }
        if ($target==3)
        {print $O3 "$line\n";
            $nb3++;
        }
        if ($target==4)
        {print $O5 "$line\n";
            $nbp++;
        }
        
	}
print "unmapped\t$sufix1\t$sufix2\tboth\tPlasmid\n";
print "$unmapped\t$nb1\t$nb2\t$nb3\t$nbp\n";
close $CMI;
close $O1;
close $O2;
close $O3;
close $O4;
close $O5;

exit;

open($O3,'<',$outfile3) or die("open: $!");
open($O1,'>>',$outfile1) or die("open: $!");
open($O2,'>>',$outfile2) or die("open: $!");

#find proba to dsitribute the reads.
my $proba=$nb1/($nb1+$nb2);

$nb1=0;
$nb2=0;
while(defined ($line = <$O3>))
{
chomp $line;
@b=split(/\t/,$line);
if (rand()<$proba)
{
  print $O1 "$b[0]\n$b[9]\n+\n$b[10]\n";
    $nb1++;
}
else
{
print $O2 "$b[0]\n$b[9]\n+\n$b[10]\n";
    $nb2++;
}
}
print "$sufix1\t$sufix2\n";
print "$nb1\t$nb2\n";

system("rm $outfile3\n");
close $O1;
close $O2;
close $O3;

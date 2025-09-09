use strict;
use warnings;
print("Requirement perl HaplotypeSample.pl SampleFile ReferencePanelFastaFile PositionsToAlinmentFile ListOfSites threads AlignmentLength");
#perl HaplotypingSample.pl ListSamples.txt MultiGenomeContig.fasta output_positions.txt SelectedSites.txt 20 251008510
if ((@ARGV) !=6)
{
    print("missing arguments\n");
    die;
}
my $penaltyScore=4;
my $matchScore=1;
my $totalpenalty=$penaltyScore+$matchScore;
my $mismatchratio=0.1;

my $SampleFile=$ARGV[0];
my $ReferencePanel=$ARGV[1];
my $PositionsToAlinmentFile=$ARGV[2];
my $ListOfSites=$ARGV[3];
my $threads=$ARGV[4];
my $Alnlen1=$ARGV[5]+1; # Added AlignmentLength as an argument

die("The max expected position must be a positive number\n")
    unless ($Alnlen1 =~ /^\d+$/ && $Alnlen1 > 0);

my $line;
my $CMI;
my $OUT;
my $SAMPLEFILES;
my (@b,$name,$Alnlen,$unmapped,$contig,$state,$alnpos,$basenumber,$genomepos);

$Alnlen=$Alnlen1;
my @Cigar;
my %CigarRead=("M"=>1,"I"=>1,"S"=>1,"D"=>0,"H"=>0);
my %CigarGenome=("M"=>1,"I"=>0,"S"=>0,"D"=>1,"H"=>0);
my %GenomeListToCluster=();
$unmapped=0;
########################possibility of new prog
#First get the list of samples


#Second read the file of position matches
open($CMI,'<',$PositionsToAlinmentFile) or die("could not the PositionsToAlinmentFile open: $!");
#open($OUT,'>',"output_positions2.txt") or die("open: $!");
my %ContigToPos=();
while(defined ($line = <$CMI>))
{
    chomp $line;
    @b=split(/\s/,$line);
    $name=$b[0];
    shift(@b);
    $ContigToPos{$name}=[@b];
    #print("@{$ContigToPos{$name}}[0..10]\n");
}
close $CMI;
#close $OUT;


#Print output either on all positons for first checks then only on choosen sites.
open($CMI,'<',$ListOfSites) or die("could not the ListOfSites ($ListOfSites) open: $!");
$line = <$CMI>; #reade header;
my $outputname="FrequenciesOutputs.txt";
open($OUT,'>',$outputname);
print $OUT "RefGenomePos\n";
while(defined ($line = <$CMI>))
    {
    chomp $line;
    @b=split('\s',$line);
    my $refGenomePos=$b[3];
    print $OUT "$refGenomePos\n";
    }
close  $CMI;
close $OUT;






open($SAMPLEFILES,'<',$SampleFile) or die("could not the SampleFile open: $!");
my $linelarge;
#Here Loop on all files
while (defined ($linelarge = <$SAMPLEFILES>))
{
chomp $linelarge;
my @c=split('\t',$linelarge);
if (@c!=3) ## This was modified as now the list of samples only has 2 columns
{
print ("missing info in $SampleFile, $linelarge\n");
die;
}
my $prefix=$c[0];

    
#system("./bwa index $ReferencePanel\n");
system("bwa index $ReferencePanel\n");

#Here loop on all samples
#Third perform Mapping against these contigs
printf("bwa mem -t $threads $ReferencePanel $c[1] $c[2] >Sample.sam\n");
system("bwa mem -t $threads $ReferencePanel $c[1] $c[2] >Sample.sam\n"); ## This was modified to use just the concatenated file
#Read the sam and output the frequency of bases at all positions

my @MatchingBases=([(0)x$Alnlen],[(0)x$Alnlen],[(0)x$Alnlen],[(0)x$Alnlen],[(0)x$Alnlen]);
my %basetonumber=("A"=>0,"C"=>1,"G"=>2, "T"=>3, "N"=>4);
my %reversebasetonumber=("A"=>3,"C"=>2,"G"=>1, "T"=>0, "N"=>4);
open($CMI,'<',"Sample.sam") or die("open: $!");
my $count=0;
while(defined ($line = <$CMI>))
{
    chomp $line;
    if ($count % 100000==0)
    {print("$count\n");}
    $count++;
    #pass headers
    if (substr($line,0,1)eq '@')
    {
        next;
    }
    #store line in array
    @b=split(/\t/,$line);
	
	
	my $readlenght=length($b[9]);
	
    $line=~/AS:i:(\d+)/;
	my $score=$1;
	
 if ($b[2] eq "*" || $score==0 )
    {
        $unmapped=$unmapped+1;
        next;
    }
	

	my $mismatch=0;

	if ($line=~/NM:i:(\d+)/)

        {$mismatch=$1;}


	my $totalM=($score+$totalpenalty*$mismatch)/$matchScore;

	my $ratio=$mismatch/$totalM;
	my $lengthratio=$totalM/$readlenght;
    #pass non mapped reads
    if ( $ratio > $mismatchratio || $lengthratio<0.8)
    {
        $unmapped=$unmapped+1;
        next;
    }
    
    #print("$b[0]\t$1\n");
    
    #extract contig name nd position in the genome
    $contig=$b[2];
    $genomepos=$b[3];
    my $readpos=0;
    #my @transit=@{$ContigToPos{$contig}};
    my $transit=$ContigToPos{$contig};
    #print("here : @transit\n");
    #get the cigar score
    @Cigar = ($b[5] =~ /(\d+\D)/g);
	#print("$contig");
    #print("$b[5]\n");
    my $kl=@{$transit};
    #print("$count\t$genomepos\t$contig\t${$transit}[$genomepos]\t$kl\n");

    foreach $state (@Cigar)
    {
        my @cigar=($state=~ /(\d+)(\D)/);
        if ($cigar[1] eq 'M')
        {
            for my $i (1..$cigar[0])
            {
                
                $alnpos=${$transit}[$genomepos];
				#print("$contig\n");
				#print("$genomepos\n");
				#print("$alnpos\n");
                if ($alnpos!=0)#base is on the alignment
                {
                    if ($alnpos>0)
                    {
                        $basenumber=$basetonumber{substr($b[9],$readpos,1)};
                    }
                    else
                    {
                         $basenumber=$reversebasetonumber{substr($b[9],$readpos,1)};
                    }
                    $alnpos=abs($alnpos);
                    $MatchingBases[$basenumber][$alnpos]=$MatchingBases[$basenumber][$alnpos]+1;
                    $genomepos++;
                    $readpos++;
                }
            }
            
        }
        else
        {
            $readpos+=$CigarRead{$cigar[1]}*$cigar[0];
			#print("$readpos\n");
            $genomepos+= $CigarGenome{$cigar[1]}*$cigar[0];
			#print("$basenumber\n");
        }
        
    }
    
}
close $CMI;

    

        
#HEre add the selection of sites for printing
my $outputname="$prefix"."_bases.txt";
open($OUT,'>',$outputname);
print $OUT "Sample\tAlnPos\tRefGenomePos\tRefAlleleCount\tAlternativeAlleleCount\tOtherCount\tRefAllele\tAlternativeAllele\tA\tC\tG\tT\tN\n";

my $inputname2="FrequenciesOutputs.txt";
open(my $IN2,'<',$inputname2) or die("could not the ListOfSites open: $!");
 
my $outputname2="FrequenciesOutputsTmp_$prefix.txt";
#my $outputname2="FrequenciesOutputsTmp.txt";
open(my $OUT2,'>',$outputname2) or die("could not the ListOfSites open: $!");
$line = <$IN2>;
chomp $line;
print $OUT2 "$line\t$prefix\n";
    
open($CMI,'<',$ListOfSites) or die("could not the ListOfSites ($ListOfSites) open: $!");
$line = <$CMI>; #reade header;

my ($pos,$allele1,$allele2,$count1 ,$count2,%seen,$count_other,$refGenomePos);
while(defined ($line = <$CMI>))
{
    chomp $line;
    @b=split('\s',$line);
    $pos=$b[0]+1;
    $allele1=$b[1];
    $allele2=$b[2];
    $refGenomePos=$b[3];
    $count1 = $MatchingBases[$basetonumber{$allele1}][$pos];
    $count2 = $MatchingBases[$basetonumber{$allele2}][$pos];
    %seen = ( $allele1 => 1, $allele2 => 1 );
    $count_other=0;
    for my $base (keys %basetonumber) {
        next if exists $seen{$base};
        $count_other += $MatchingBases[$basetonumber{$base} ][$pos];
    }
    print $OUT   "$prefix\t$pos\t$refGenomePos\t$count1\t$count2\t$count_other\t$allele1\t$allele2\t".$MatchingBases[0][$pos]."\t".$MatchingBases[1][$pos]."\t".$MatchingBases[2][$pos]."\t".$MatchingBases[3][$pos]."\t".$MatchingBases[4][$pos]."\n";
    
    $line = <$IN2>;
    chomp $line;
    my $sum=$count1+$count2+$count_other;
   
    my $formatted=0;
    if ($sum>0)
    {my $freq=$count2/$sum;
     $formatted = sprintf("%.3f", $freq);
    }else{$formatted ="NA";}
    
    print $OUT2 "$line\t$formatted\n";
}
close  $CMI;
close $OUT;
close $OUT2;
close $IN2;
system("rm FrequenciesOutputs.txt\n");
system("mv FrequenciesOutputsTmp.txt FrequenciesOutputs.txt\n");
#    while(defined ($line = <$CMI>))
#{
#chomp $line;
#@b=split(',',$line);

#my $outputname="$prefix"."_$b[0]"."_bases.txt";
#open($OUT,'>',$outputname);
#if ($b[0] eq "all")
#{
#@b=(1..($Alnlen-1))
##}
#else
#{
#shift(@b);
#}
#for my $i (@b)
#{
#print $OUT "$prefix\t$i\t".$MatchingBases[0][$i]."\t".$MatchingBases[1][$i]."\t".$MatchingBases[2][$i]."\t".$MatchingBases[3][$i]."\t".$MatchingBases[4][$i]."\n";
#}
#close $OUT;
#}
system ("rm Sample.sam");
}
close $SAMPLEFILES;
die;


#store aln and filter positions (keep only variable sites)

#this code takes an xmfa file form parnsp alginelmtn from two genomes and exclude the blocs wher efor the second genome the distance to the previous block is more than 100kb and exclude also all the small blocks of less than 110 bp.
#perl CreateFilesForHaplotypingSamples_Recombination.pl parsnp.xmfa ../../../Ref_Genomes/536
#!/usr/bin/perl
use strict;
use warnings;
my $CMI;
#my $input = 'parsnp.xmfa';
# Check arguments
die "Usage: $0 <parsnp.xmfa> <genome_directory>\n" unless @ARGV == 2;
#my ($xmfa) = @ARGV;
my ($xmfa, $genome_dir) = @ARGV;
my $maxdist=100000;
my $minlen=110;

#get the fasta fils names from the xmfa file
my (%files, %declared_len);

# Step 1: Parse header
open my $fh, '<', $xmfa or die "Cannot open $xmfa: $!";
while (my $line = <$fh>) {
    chomp $line;
    last if $line =~ /^#IntervalCount/;
    if ($line =~ /^##SequenceIndex\s+(\d+)/) {
        my $index = $1;
        while (my $subline = <$fh>) {
            chomp $subline;
            last if $subline =~ /^##SequenceIndex/;
            last if $subline =~ /^#IntervalCount/;
            if ($subline =~ /^##SequenceFile\s+(.+)/) {
                my $file = $1;
                $file =~ s/\s+//g;
                $files{$index} = $file;
            }
            if ($subline =~ /^##SequenceLength\s+(\d+)bp/) {
                $declared_len{$index} = $1;
                last;  # fini la lecture de ce bloc SequenceIndex
            }
        }
    }
}
close $fh;

close $fh;
use Data::Dumper;
print Dumper(\%files);

#get the gnome length
sub get_fasta_length {
    my ($file) = @_;
    open my $f, '<', $file or die "Cannot open $file: $!";
    my $len = 0;
    my $seq = '';
    while (<$f>) {
        chomp;
        next if /^>/;
        $seq .= $_;
        $len += length($_);
    }
    close $f;
    return ($len, $seq);
}

# Step 3: Report
open($CMI,'>',"MultiFastaforBWA.fasta") or die("open: $!");

foreach my $index (sort { $a <=> $b } keys %files) {
    my $filename = $files{$index};
    my $fullpath = "$genome_dir/$filename";
    my ($computed_len,$fastagenomeseq) = get_fasta_length($fullpath);
    print $CMI ">$files{$index}\n$fastagenomeseq\n";
#    print "Genome $index: $filename\n";
#   print " - Declared length : $declared_len{$index} bp\n";
 #   print " - Computed length : $computed_len bp\n";
}
close $CMI;

open $fh, '<', $xmfa or die "Cannot open $xmfa: $!";

my @blocks;
my @current;

while (my $line = <$fh>) {
    chomp $line;

    if ($line =~ /^>/) {
        # nouveau bloc : on stocke le bloc précédent s'il existe
        push @blocks, [@current] if @current;
        @current = ([$line, '']);   # ligne d'entête + séquence vide
    } elsif (@current && defined $current[-1][1]) {
        # ajouter la séquence au dernier bloc courant
        $current[-1][1] .= $line;
    }
}


# ajouter le dernier bloc
push @blocks, [@current] if @current;

print "Nombre de blocs totaux : ", scalar(@blocks), "\n";


# Parse and rank blocks by genome1 start position
my @parsed_blocks;

for my $block (@blocks) {
    next unless @$block == 2;

    my @entries;
    foreach my $entry (@$block) {
        my ($header, $seq) = @$entry;
        my ($seqid, $start, $end, $strand) = $header  =~/^>(\d+):(\d+)-(\d+) ([+-])/;
        push @entries, { seqid=>$seqid, start => $start + 0, end => $end + 0, strand => $strand, seq => $seq };
    }

    # Sort entries by start → genome1 assumed to be first
    my ($g1, $g2) = @entries;#sort { $a->{start} <=> $b->{start} } @entries;

    push @parsed_blocks, {
        g1_id=> $g1->{seqid},
        g2_id=> $g2->{seqid},
        g1_start => $g1->{start},
        g2_start => $g2->{start},
        g1_end => $g1->{end},
        g2_end => $g2->{end},
        g2_strand => $g2->{strand},
        seq1 => $g1->{seq},
        seq2 => $g2->{seq}
    };
}

# Sort blocks by genome1 start
@parsed_blocks = sort { $a->{g1_start} <=> $b->{g1_start} } @parsed_blocks;

# Iteratively select syntenic blocks
my @syntenic;
my ($prev_end, $prev_strand);

for my $i (0 .. $#parsed_blocks) {
    my $blk = $parsed_blocks[$i];
    if ($i == 0) {
        push @syntenic, $blk;
        $prev_end  = $blk->{g2_end};
        $prev_strand = $blk->{g2_strand};
        next;
    }

    my $delta = abs($blk->{g2_start} - $prev_end);
    my $len=$blk->{g2_end}-$blk->{g2_start};
    if ($blk->{g2_strand} eq $prev_strand && $delta <= $maxdist && $len>$minlen) {
        push @syntenic, $blk;
        
        $prev_end  = $blk->{g2_end};
        $prev_strand = $blk->{g2_strand};
        
        #print("$blk->{g1_start}\t$blk->{g2_start}\t$len\t$delta\n");
    }
    # else: skip non-syntenic block
}

# Reconstruct alignment per genome
my ($genome1, $genome2) = ('', '');
#my $genome_len = 5000000;  # or however long your genome is
my @genomecorrespondance1 = (0) x ($declared_len{1}+1);
my @genomecorrespondance2 = (0) x ($declared_len{2}+1);


my @syntenyruptures;
my $AlignmentLength=0;
my $aln_pos=0;
for my $blk (@syntenic) {
    my $alnposstart=$aln_pos;
    my $aln_pos_mem=$aln_pos;
    #ref strand is alway on plus strand
    my $seq   = $blk->{seq1};
    my $g1_index = $blk->{g1_start};  # genome coordinate (0-based or 1-based!)
    #printf "$blk->{g1_id}\t$g1_index\n";
    for (my $i = 0; $i < length($seq); $i++) {
        my $base = substr($seq, $i, 1);
            if ($base ne '-') {
                $genomecorrespondance1[$g1_index] = $aln_pos;
                $g1_index++;
            }
            $aln_pos++;
        }
    
    #rest $alnpos
    $aln_pos=$aln_pos_mem;
    $seq      = $blk->{seq2};
    if ($blk->{g2_strand} eq "+")
    {
        my $g2_index = $blk->{g2_start};  # genome coordinate (0-based or 1-based!)
        for (my $i = 0; $i < length($seq); $i++) {
            my $base = substr($seq, $i, 1);
            if ($base ne '-') {
                $genomecorrespondance2[$g2_index] = $aln_pos;
                $g2_index++;
            }
            $aln_pos++;
        }
    }else
    {
        my $g2_index = $blk->{g2_end};  # genome coordinate (0-based or 1-based!)
        for (my $i = 0; $i < length($seq); $i++) {
            my $base = substr($seq, $i, 1);
            if ($base ne '-') {
                $genomecorrespondance2[$g2_index] = -$aln_pos;
                $g2_index--;
            }
            $aln_pos++;
        }
    }
    
    
    $genome1 .= $blk->{seq1};
    $genome2 .= $blk->{seq2};
    #soquer les points de rupture....
    my $alnposend=$aln_pos-1;
    push @syntenyruptures, {start=>$alnposstart,end=>$alnposend,refstart=>$blk->{g1_start}, refend=>$blk->{g1_end}} ;
}
$AlignmentLength=$aln_pos;
print "ALn Length: $AlignmentLength\n";


#Now the reverse find in the reference genome the adreess of the bases found in the aln.
$aln_pos=0;
my @alncorrespondance = (0) x ($AlignmentLength+1);
for my $blk (@syntenic) {
    my $seq   = $blk->{seq1};
    my $g1_index = $blk->{g1_start};  # genome coordinate (0-based or 1-based!)
    #printf "$blk->{g1_id}\t$g1_index\n";
    for (my $i = 0; $i < length($seq); $i++) {
        my $base = substr($seq, $i, 1);
        $alncorrespondance[$aln_pos] = $g1_index;
        if ($base ne '-') {
            #$alncorrespondance[$aln_pos] = $g1_index;
            $g1_index++;
        }
        $aln_pos++;
    }
    }


open( $CMI,'>',"output_positions.txt") or die("open: $!");
print $CMI "$files{1}"." @genomecorrespondance1\n";
print $CMI "$files{2}"." @genomecorrespondance2\n";
close $CMI;

#now find the polymrophic sites
my $len = length($genome1);
my $genomerefpos;
die "Aligned sequences not same length\n" unless $len == length($genome2);
open( $CMI,'>',"AlnPolymorphismPositions.fasta") or die("open: $!");
print $CMI "Position_in_Aln\t$files{1}\t$files{2}\tPosition_in_ref_file\n"; #TO BE PRINTED IN A FILE

my $blkindex=0;
my $start = $syntenyruptures[$blkindex]{start};
my $end   = $syntenyruptures[$blkindex]{end};
for (my $i = 0; $i < $len; $i++) {
    if ($i>$syntenyruptures[$blkindex]{end})
    {$blkindex=$blkindex+1;
    $start = $syntenyruptures[$blkindex]{start};
    $end   = $syntenyruptures[$blkindex]{end};
    }
    
    my $a1 = substr($genome1, $i, 1);
    my $a2 = substr($genome2, $i, 1);

    next if $a1 eq $a2;
    #
    #store print the distance to eitehr previous polym, or previous rupture
    splice(@syntenyruptures, $blkindex, 1);
    my @new_elements;
    my $rs=$alncorrespondance[$start];
    my $re=$alncorrespondance[$i-1];
    if ($rs*$re==0)
    {print ("$start\t$rs\t$i\t$re\t$a1\t$a2\n");}
    if ($i > $start)
        {
        push @new_elements, { start => $start, end => $i - 1 , refstart=>$alncorrespondance[$start], refend=>$alncorrespondance[$i-1]};
        }
    if ($i<$end)
        {
        push @new_elements, { start => $i + 1, end => $end   , refstart=>$alncorrespondance[$i+1], refend=>$alncorrespondance[$end] } ;
        }
    
    splice(@syntenyruptures, $blkindex, 0, @new_elements);
    $start = $syntenyruptures[$blkindex]{start};
    $end   = $syntenyruptures[$blkindex]{end};
    
    
    next if $a1 eq '-' || $a2 eq '-';
    $genomerefpos=$alncorrespondance[$i];
    print $CMI "$i\t$a1\t$a2\t$genomerefpos\n"; #TO BE PRINTED IN A FILE
}
close $CMI;
# Output FASTA
open( $CMI,'>',"FilteredAln.fasta") or die("open: $!");
print $CMI ">$files{1}\n$genome1\n";#TO BE PRINTED IN A FILE
print $CMI ">$files{2}\n$genome2\n";
close $CMI;

print "Nombre de blocs totaux : ", scalar(@parsed_blocks), "\n";
print "Nombre de blocs synteniques retenus : ", scalar(@syntenic), "\n";


open  $fh, '>', 'PerfectHomologyBlocs.txt' or die "Cannot open file: $!";
print $fh "Aln_start\tAln_end\tRef_start\tRef_end\n";
foreach my $rupture (@syntenyruptures) {
    print $fh "$rupture->{start}\t$rupture->{end}\t$rupture->{refstart}\t$rupture->{refend}\n";
}
close $fh;

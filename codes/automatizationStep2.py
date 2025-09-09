import os
import subprocess

#os.environ["PATH"] = "/shared/ifbstor1/software/miniconda/envs/bwa-0.7.17/bin:" + os.environ.get("PATH", "")


folder = "data/qualityreduced175"
perl_file = "PaperConj2025_files/HaplotypingSampleRecombination_040925.pl"
SampleFile =  "PaperConj2025_files/Step1/536_0175/ListSamples.txt"
ReferencePanelFastaFile = "PaperConj2025_files/parsnp536/MultiFastaforBWA.fasta"
PositionsToAlignmentFile = "PaperConj2025_files/parsnp536/output_positions.txt" 
ListOfSites = "PaperConj2025_files/parsnp536/AlnPolymorphismPositions.fasta"
nbthreads = 12
AlignmentLength = 3939501  


subprocess.run(["perl", perl_file, 
                        SampleFile, 
                        ReferencePanelFastaFile, 
                        PositionsToAlignmentFile,
                        ListOfSites, 
                        str(nbthreads), 
                        str(AlignmentLength)])
                        
#perl "PaperConj2025_files/HaplotypingSampleRecombination2.pl" "PaperConj2025_files/Step1/ListSamples.txt" "PaperConj2025_files/Step1/MultiFastaforBWA.fasta" "PaperConj2025_files/Step1/output_positions.txt" "PaperConj2025_files/Step1/AlnPolymorphismPositions.fasta" 20 251008510

#changement du dernier index pas le txt mais le fasta Aln 
## 👩‍💻 Codes

This folder contains the following files:  

📂 codes/  
├── 📄 PaperConj2025_analysis_090925.Rmd  
├── 📄 bwa_25.py  
├── 📄 analysis_downscaling.py  
├── 📄 downscaling.sh  
├── 📄 ecoli_25.py  
└── 📄 SortReads.pl  

____

**R files**
- `PaperConj2025_analysis_090925.Rmd`    
This .R notebook regroup all the codes needed for the paper main and supplementary figures. The file uses different inputs that are either the supplementary tables or output files (e;g. .txt, .sam..) of previous perl/python pipeline.


**Bash files**
- `downscaling.sh`    
Using seqtk lower the sequencing resolution of chosen fastq files. 


**python files**
- `bwa_25.py`
- 'analysis_downscaling.py' calculated the coverage of fastq files
- `ecoli_25.py`  
<small> These python files are used for .fastaq aligned by BWA using `bwa_25.py` and their parsing automatized in a `ecoli_25.py` </small>

**perl files**  
- `SortReads.pl`  
<small> This .pl is an alternative to `ecoli_25.py` after BWA aligment. </small>

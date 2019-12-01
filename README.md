# CAMAMED
A pipeline for composition-aware mapping-based analysis of metagenomic data
Software Requirements
-	Linux operating system (Preferably Ubuntu)
-	Python 2.7
-	MetaPhlAn2 
    -Installation command
        -sudo apt install metaphlan2 
    -After the first run, MetaPhlAn database files are downloaded automatically. Otherwise, download the files from the following links        and copy them to the installation path in folder /usr/share/metaphlan2/databases. 
        •	https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar 
        •	https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.md5
    -	CD-HIT
       -Installation command
           -sudo apt-get install cd-hit
    -	SRA-Toolkit
      -Installation command
          -sudo apt install sra-toolkit
    -	Samtools
      -Installation command
          -sudo apt-get install samtools
    -	Also install the necessary packages for Python if requested. For example:
          -sudo apt-get install python-numpy
          -sudo apt-get install python-pandas
          -pip install biopython
    -	Install R package for python 2.7 
          -pip install rpy2==2.8.6
          -conda install -c r r-glmnet 
          -conda install -c r r-kernsmooth
          -python2
          -from rpy2.robjects.packages import importr
          -base = importr('base')
          -base.source("http://bioconductor.org/biocLite.R")
          -biocinstaller = importr("BiocInstaller")
          -biocinstaller.biocLite("metagenomeSeq")
          -biocinstaller.biocLite("edgeR")

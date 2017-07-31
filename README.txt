------------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Web server: http://apoplastp.csiro.au/

Software download and installation: http://apoplastp.csiro.au/software.html

Instructions: http://apoplastp.csiro.au/instructions.html
------------------------------------------------------------------------------
------------------------------------------------------------------------------
------------------------------------------------------------------------------

------------------------------------------------------------------------------
------------------------------------------------------------------------------
ApoplastP: prediction of effectors and plant proteins in the apoplast using machine learning
Copyright (C) 2017-2018 Jana Sperschneider	
Contact: jana.sperschneider@csiro.au
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Installation instructions for ApoplastP 1.0
------------------------------------------------------------------------------
------------------------------------------------------------------------------
ApoplastP relies on two tools, the EMBOSS software and the WEKA 3.8.1 software. These have been shipped 
with the ApoplastP 1.0.tar.gz file, but they need to be installed by the user. 
-----------------------------------------
1) Make sure ApoplastP has the permission to execute. Then unpack ApoplastP in your desired location:

chmod -R 755 ApoplastP_1.0.tar.gz
tar xvf ApoplastP_1.0.tar.gz
cd ApoplastP_1.0 
-----------------------------------------

2) Install EMBOSS

-----------------------------------------
cd Scripts
tar xvf emboss-latest.tar.gz 
cd EMBOSS-6.5.7/
./configure
make
cd ../
-----------------------------------------

3) Install WEKA: simply unzip the file weka-3-8-1.zip

-----------------------------------------
unzip weka-3-8-1.zip
-----------------------------------------

3) Run ApoplastP

To test that ApoplastP is working, type the following command in the working directory ApoplastP_1.0/Scripts

-----------------------------------------
python ApoplastP.py -i Testing.fasta
-----------------------------------------

 If you are getting an error message like 'ImportError: No module named Bio', you need to install BioPython on your computer. See here for help. For example, you can try and run 'pip install biopython'. 
Note that ApoplastP runs under Python 2.x, not under Python 3.x.
If you are having troube installing EMBOSS, please see here for help: http://emboss.sourceforge.net/download/
If you are having troube installing WEKA, please see here for help: http://www.cs.waikato.ac.nz/~ml/weka/index.html


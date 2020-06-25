#### What is ApoplastP?

The plant apoplast is vital to signalling, transport and plant-pathogen interactions. 
 ApoplastP is a machine learning method for predicting localization of proteins to the plant apoplast.
 ApoplastP can distinguish non-apoplastic proteins from apoplastic proteins for both plant proteins and pathogen proteins.
 In particular, ApoplastP can predict if an effector localizes to the plant apoplast. ApoplastP achieves sensitivity of around 75% and specificity of over 95% on test sets of effectors. 
 
#### What is ApoplastP not?

ApoplastP is not a tool for secretome prediction.

ApoplastP has been trained to predict proteins that localize to the plant apoplast, so please submit a FASTA file of secreted proteins to test if they are predicted apoplastic proteins. It is recommended to use tools such as SignalP or Phobius	to predict first if a protein is likely to be secreted. Alternatively, experimentally determined secretomes (e.g. apoplastic proteomics) instead of computationally predicted secretomes can be submitted to ApoplastP to discover unconventionally secreted proteins that localize to the apoplast. 

#### Installing ApoplastP

You can submit your proteins of interest to the webserver at http://apoplastp.csiro.au/ or install it locally.
All training and evaluation data can be found [here](http://apoplastp.csiro.au/data.html).

ApoplastP has been written in Python and uses pepstats from the EMBOSS software and the WEKA 3.6 software. It also requires that you have BioPython installed. **ApoplastP from version 1.0.2 inclusive uses Python 3.** 

To get ApoplastP to work on your local machine, you need to install the EMBOSS and WEKA 3.8.1 softwares from source. Both are already provided in the ApoplastP distribution to ensure that compatible versions are used. 

0. Download the latest release from this github repo (or alternatively you can clone the github repo and skip step 1).

1. Make sure ApoplastP has the permission to execute. Then unpack LOCALIZER in your desired location
```
tar xvf ApoplastP_1.0.2.tar.gz
chmod -R 755 ApoplastP_1.0.2/
cd ApoplastP_1.0.2
```

2. For the EMBOSS installation, you need to switch to the Scripts directory and unpack, configure and make. Alternatively, if you are on a computer cluster and EMBOSS is already installed, you can change the variable PEPSTATS_PATH in the ApoplastP.py script to the EMBOSS directory that contains pepstats on the machine you are using.
```
cd Scripts
tar xvf emboss-latest.tar.gz
cd EMBOSS-6.5.7/
./configure
make
cd ../ 
```

3. For WEKA, you need to simply unzip the file unzip weka-3-8-1.zip
```
unzip weka-3-8-1.zip
```
If you are having troube installing EMBOSS, please see [here](http://emboss.sourceforge.net/download/) for help.
If you are having troube installing WEKA, please see [here](https://www.cs.waikato.ac.nz/~ml/weka/index.html) for help. 

4. Test if ApoplastP is working
```
python ApoplastP.py -i Testing.fasta
```

5. Problems?

If you are getting an error message like 'ImportError: No module named Bio', you need to install BioPython on your computer. See [here](https://biopython.org/wiki/Download) for help. For example, you can try and run:
```
pip install biopython
```

#### ApoplastP output format

Run this to get a feel for the output format:
```
python ApoplastP.py -i Testing.fasta

-----------------

ApoplastP is running for 29 proteins given in FASTA file Testing.fasta

-----------------
Calculate statistics of protein properties
# Identifier     Prediction      Probability
Ecp6    Apoplastic      0.96
AvrLm6  Apoplastic      0.78
AvrLm1  Non-apoplastic  0.51
AvrLm11 Apoplastic      0.56
Nep1    Apoplastic      0.76
ToxB    Apoplastic      0.65
Msp1    Apoplastic      0.87
AvrStb6 Apoplastic      0.76
Cgfl    Apoplastic      0.81
PstSCR1 Non-apoplastic  0.54
CfTom1  Non-apoplastic  0.69
AvrLm4-7        Apoplastic      0.82
Ave1    Apoplastic      0.79
AVR9    Apoplastic      0.8
AVR4    Apoplastic      0.93
AVR4E   Apoplastic      0.72
AVR2    Apoplastic      0.83
Avr5    Apoplastic      0.92
Ecp2    Apoplastic      0.87
Ecp1    Apoplastic      0.96
Ecp5    Apoplastic      0.88
Ecp4    Apoplastic      0.96
Bas4    Apoplastic      0.82
MC69    Apoplastic      0.85
Slp1    Apoplastic      0.98
NIP1    Apoplastic      0.91
Pep1    Apoplastic      0.95
Pit2    Apoplastic      0.85
Mg3LysM Apoplastic      1.0

-----------------
Predicted apoplastic proteins:

Ecp6| Apoplastic probability:0.96
AvrLm6| Apoplastic probability:0.78
AvrLm11| Apoplastic probability:0.56
Nep1| Apoplastic probability:0.76
ToxB| Apoplastic probability:0.65
Msp1| Apoplastic probability:0.87
AvrStb6| Apoplastic probability:0.76
Cgfl| Apoplastic probability:0.81
AvrLm4-7| Apoplastic probability:0.82
Ave1| Apoplastic probability:0.79
AVR9| Apoplastic probability:0.8
AVR4| Apoplastic probability:0.93
AVR4E| Apoplastic probability:0.72
AVR2| Apoplastic probability:0.83
Avr5| Apoplastic probability:0.92
Ecp2| Apoplastic probability:0.87
Ecp1| Apoplastic probability:0.96
Ecp5| Apoplastic probability:0.88
Ecp4| Apoplastic probability:0.96
Bas4| Apoplastic probability:0.82
MC69| Apoplastic probability:0.85
Slp1| Apoplastic probability:0.98
NIP1| Apoplastic probability:0.91
Pep1| Apoplastic probability:0.95
Pit2| Apoplastic probability:0.85
Mg3LysM| Apoplastic probability:1.0
-----------------

Number of proteins that were tested: 29
Number of predicted apoplastic proteins: 26

-----------------
89.7 percent are predicted to be apoplastic.
-----------------
```

ApoplastP will return the output as shown in the example above. First, it will return the predicted apoplastic proteins in your set as FASTA sequences, if there are any. Second, a summary table will be shown which shows the predictions (apoplastic or non-apoplastic) for each submitted protein.

ApoplastP returns a probability that a tested instance will belong to either the apoplastic or non-apoplastic class.

We deliberately did not recommend a probability threshold over which a protein would be classified as an apoplastic protein candidate, as we believe it should remain up to the individual user to interpret their results in the context of additional resources available. For example, a researcher might like to predict the full apoplastic candidate list and overlay this with in planta expression data to prioritize candidates, whereas in other situations without additional information a list of high-priority candidates as determined by the ApoplastP probabilities might be more appropriate. 

#### Citation for ApoplastP:

Sperschneider J et al. (2017) ApoplastP: prediction of effectors and plant proteins in the apoplast using machine learning. New Phytologist. doi:10.1111/nph.14946. 

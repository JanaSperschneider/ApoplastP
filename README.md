#### What is ApoplastP?

The plant apoplast is vital to signalling, transport and plant-pathogen interactions. 
 ApoplastP is a machine learning method for predicting localization of proteins to the plant apoplast.
 ApoplastP can distinguish non-apoplastic proteins from apoplastic proteins for both plant proteins and pathogen proteins.
 In particular, ApoplastP can predict if an effector localizes to the plant apoplast. ApoplastP achieves sensitivity of around 75% and specificity of over 95% on test sets of effectors. 
 
#### What is ApoplastP not?

ApoplastP is not a tool for secretome prediction.

ApoplastP has been trained to predict proteins that localize to the plant apoplast, so please submit a FASTA file of secreted proteins to test if they are predicted apoplastic proteins. It is recommended to use tools such as SignalP or Phobius	to predict first if a protein is likely to be secreted. Alternatively, experimentally determined secretomes (e.g. apoplastic proteomics) instead of computationally predicted secretomes can be submitted to ApoplastP to discover unconventionally secreted proteins that localize to the apoplast. 

#### Citation for ApoplastP:

Sperschneider J et al. (2017) ApoplastP: prediction of effectors and plant proteins in the apoplast using machine learning. New Phytologist. doi:10.1111/nph.14946. 

#### Running ApoplastP

You can submit proteins to the webserver at http://apoplastp.csiro.au/.

Alternatively, you can install ApoplastP on your machine to run it locally.
For detailed installation instructions see here: http://apoplastp.csiro.au/software.html

For help on how to interpret the output format, see http://apoplastp.csiro.au/output.html

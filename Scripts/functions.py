import os
import sys
import re
import subprocess

import os
import sys
import re
import random
import numpy as np
#from scipy.spatial.distance import pdist, squareform
import subprocess as sub

from Bio import SeqIO
from Bio.SeqUtils import ProtParam
import os
import sys
import re
import io
import getopt
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()

        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running EffectorP.       
    """
    print '''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# EffectorP :: predicting fungal effector proteins from secretomes using machine learning
# EffectorP 1.0 (July 2015); http://effectorp.csiro.au/
# Copyright (C) 2015-2016 Jana Sperschneider, CSIRO.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''
    print "Usage for EffectorP: ", 
    print "python EffectorP.py [-options] -i <input_file>"
    print 
    print "where basic options are:"
    print "-h : show brief help on version and usage" 
    print 
    print "options for output format:"
    print "-s : short output format that provides predictions for all proteins as one tab-delimited table [default long format]"
    print
    print "options directing output:"
    print "-o <f> : direct output to file <f>, not stdout"
    print "-E <f> : save predicted effectors to FASTA file <f>"        
    print
    print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
    print
    sys.exit()    

    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()

        Purpose:  Scan the input options given to the EffectorP program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "hso:E:i:", ["help"])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    FASTA_FILE = None
    short_format = False
    output_file = None
    effector_output = None

    i_count, o_count, E_count, P_count = 0, 0, 0, 0
   
    for opt, arg in opts:
        if opt in ("-o"):
            output_file = arg
            o_count += 1
        elif opt in ("-s"):
            short_format = True
        elif opt in ("-i"):
            FASTA_FILE = arg
            i_count += 1
        elif opt in ("-E"):
            effector_output = arg
            E_count += 1
        elif opt in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if i_count > 1 or o_count > 1 or E_count > 1:
       usage()

    return FASTA_FILE, short_format, output_file, effector_output
# -----------------------------------------------------------------------------------------------------------
# Taken from http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
GRAVY_DIC = {
'A': 1.8, 
'R': -4.5,
'N': -3.5,
'D': -3.5,  
'C':  2.5,  
'Q': -3.5,  
'E': -3.5,  
'G': -0.4,  
'H': -3.2,  
'I':  4.5,  
'L':  3.8,  
'K': -3.9,  
'M':  1.9,  
'F':  2.8,  
'P': -1.6,  
'S': -0.8,  
'T': -0.7,  
'W': -0.9,  
'Y': -1.3,  
'V':  4.2}
# -----------------------------------------------------------------------------------------------------------
def GRAVY(sequence):
    # The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values of all the amino acids, 
    # divided by the number of residues in the sequence. 

    gravy = 0.0
    for aa in sequence:
	    if aa.upper() in GRAVY_DIC:
	        gravy += GRAVY_DIC[aa.upper()]

    gravy = gravy/len(sequence)

    return gravy
# -----------------------------------------------------------------------------------------------------------
ARFF_APOPLAST_HEADER = '''@RELATION effectors
@ATTRIBUTE Charge NUMERIC
@ATTRIBUTE Isoelectric NUMERIC
@ATTRIBUTE Tiny NUMERIC
@ATTRIBUTE Small NUMERIC
@ATTRIBUTE Aliphatic NUMERIC
@ATTRIBUTE Aromatic NUMERIC
@ATTRIBUTE Nonpolar NUMERIC
@ATTRIBUTE Polar NUMERIC
@ATTRIBUTE Charged NUMERIC
@ATTRIBUTE Basic NUMERIC
@ATTRIBUTE Acidic NUMERIC
@ATTRIBUTE A NUMERIC
@ATTRIBUTE C NUMERIC
@ATTRIBUTE D NUMERIC
@ATTRIBUTE E NUMERIC
@ATTRIBUTE F NUMERIC
@ATTRIBUTE G NUMERIC
@ATTRIBUTE H NUMERIC
@ATTRIBUTE I NUMERIC
@ATTRIBUTE K NUMERIC
@ATTRIBUTE L NUMERIC
@ATTRIBUTE M NUMERIC
@ATTRIBUTE N NUMERIC
@ATTRIBUTE P NUMERIC
@ATTRIBUTE Q NUMERIC
@ATTRIBUTE R NUMERIC
@ATTRIBUTE S NUMERIC
@ATTRIBUTE T NUMERIC
@ATTRIBUTE V NUMERIC
@ATTRIBUTE W NUMERIC
@ATTRIBUTE Y NUMERIC
@ATTRIBUTE Gravy NUMERIC
@ATTRIBUTE Aromaticity NUMERIC
@ATTRIBUTE Instability NUMERIC
@ATTRIBUTE CysteineCount NUMERIC
@ATTRIBUTE class {APOPLAST,NON_APOPLAST}
@DATA
'''
#@ATTRIBUTE MW NUMERIC
# -----------------------------------------------------------------------------------------------------------
def write_FASTA(f_output, IDENTIFIERS, SEQUENCES):

    f = open(f_output, 'w')

    for ident, seq in zip(IDENTIFIERS, SEQUENCES):
        f.writelines(ident)
        f.writelines(seq + '\n')
    f.close()

    return
# -----------------------------------------------------------------------------------------------------------
def featurevector_apoplast(DATASET_POSITIVE, DATASET_NEGATIVE, pepstats_dic_positive, pepstats_dic_negative):

    X = []

    for TARGET_ID, sequence in DATASET_POSITIVE:
        TARGET_ID = TARGET_ID.replace('>','')
        TARGET_ID = TARGET_ID.strip()        
                    
        molecular_weight, charge, isoelectric, tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic, amino_acid_frequencies, length, dayhoff = pepstats_dic_positive[TARGET_ID]

        prot = ProtParam.ProteinAnalysis(sequence.replace('*',''))                

        #feature_vector = [charge, isoelectric] + [tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic] + amino_acid_frequencies + [GRAVY(sequence)] + [prot.aromaticity(), prot.instability_index(), prot.molecular_weight()] + [sequence.count('C')] 

        feature_vector = [charge, isoelectric] + [tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic] + amino_acid_frequencies + [GRAVY(sequence)] + [prot.aromaticity(), prot.instability_index()] + [sequence.count('C')] 

        X.append(feature_vector)

    print '-------------'
    for TARGET_ID, sequence in DATASET_NEGATIVE:
        TARGET_ID = TARGET_ID.replace('>','')
        TARGET_ID = TARGET_ID.strip()        

        molecular_weight, charge, isoelectric, tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic, amino_acid_frequencies, length, dayhoff = pepstats_dic_negative[TARGET_ID]

        prot = ProtParam.ProteinAnalysis(sequence.replace('*',''))            

        #feature_vector = [charge, isoelectric] + [tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic] + amino_acid_frequencies + [GRAVY(sequence)] + [prot.aromaticity(), prot.instability_index(), prot.molecular_weight()] + [sequence.count('C')] 

        feature_vector = [charge, isoelectric] + [tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic] + amino_acid_frequencies + [GRAVY(sequence)] + [prot.aromaticity(), prot.instability_index()] + [sequence.count('C')] 

        X.append(feature_vector)

    X = np.array(X)

    return X
# -----------------------------------------------------------------------------------------------------------
def evaluate_performance(FILENAME):

    f = open(FILENAME, 'r')
    content = f.readlines()    
    f.close()

    #TP Rate  FP Rate  Precision  Recall   F-Measure  MCC      ROC Area  PRC Area  Class
    #0.712    0.058    0.805      0.712    0.756      0.683    0.940     0.845     APOPLAST
    sensitivity, specificity, kappa, Fmeasure, MCC, AUC = None, None, None, None, None, None

    for index, line in enumerate(content):
        if '=== Stratified cross-validation ===' in line:
            for followingindex, followingline in enumerate(content[index:]):
   	        if 'Kappa statistic' in followingline:
                    kappa = float(followingline.split('Kappa statistic')[1])
                if '=== Detailed Accuracy By Class ===' in followingline:     
                    line1 = content[index + followingindex+3]
                    line2 = content[index + followingindex+4]             
                    sensitivity = float(line1.split()[0])
                    specificity = float(line2.split()[0])
                    Fmeasure = float(line1.split()[4])            
                    MCC = float(line1.split()[5])                    
                    AUC = float(line1.split()[6])                    
                    
    return sensitivity, specificity, kappa, Fmeasure, MCC, AUC
# -----------------------------------------------------------------------------------------------------------
def parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES):
    """ Function: parse_weka_output()

        Purpose:  Given the WEKA output file and the query identifiers and sequences, 
                  parse the predicted class for each protein from the WEKA output. 
                  Write the predicted effectors to a FASTA file.
              
        Input:    WEKA output file and the query identifiers and sequences.                  
    
        Return:   The set of predicted effectors only as well as all predictions. 
    """    
    predicted_apoplast, not_predicted_apoplast, predictions = [], [], []

    with open(file_input) as f:

        content = f.readlines()        

        #content_start = content.index('    inst#     actual  predicted error prediction (Charge,Isoelectric,Tiny,Small,Aliphatic,Aromatic,Nonpolar,Polar,Charged,Basic,Acidic,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Gravy,Aromaticity,Instability,MW,CysteineCount)\n')   

        content_start = content.index('    inst#     actual  predicted error prediction (Charge,Isoelectric,Tiny,Small,Aliphatic,Aromatic,Nonpolar,Polar,Charged,Basic,Acidic,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Gravy,Aromaticity,Instability,CysteineCount)\n')  
        content = content[content_start + 1:]

        for line in content:
            if line.strip():
                position = line.split()[0]
                prediction = line.split()[2]
                prob = float(line.split()[3])        
 
                # WEKA output counts from position 1, our identifiers are counted from zero
                identifier = ORIGINAL_IDENTIFIERS[int(position) - 1]
                sequence = SEQUENCES[int(position) - 1]

                if 'NON_APOP' in prediction:                               
                    noneffector = identifier.strip()
                    noneffector = noneffector.replace('>', '')  
                    predictions.append((noneffector, 'Non-apoplastic', prob, sequence))
                    not_predicted_apoplast.append((noneffector, prob, sequence))
                else:                    
                    effector = identifier.strip()
                    effector = effector.replace('>', '')            
                    #predictions.append((effector, 'Apoplastic', prob, sequence))
                    #predicted_apoplast.append((effector, prob, sequence))
                    if float(prob) >= 0.55:                                   
                        predictions.append((effector, 'Apoplastic', prob, sequence))
                        predicted_apoplast.append((effector, prob, sequence))
                    else:
                        predictions.append((effector, 'Weakly apoplastic', prob, sequence))
                        not_predicted_apoplast.append((effector, prob, sequence))

    return predicted_apoplast, not_predicted_apoplast, predictions
# -----------------------------------------------------------------------------------------------------------
def pepstats_aas(DATASET, TMP_PATH):
 
    f = open(TMP_PATH + 'temp_pepstats.fasta', 'w')

    for ident, seq in DATASET:
        f.writelines('>' + ident.replace('>','') + '\n')
        f.writelines(seq + '\n')
    f.close()

    os.popen('pepstats -sequence ' + TMP_PATH + 'temp_pepstats.fasta -outfile ' + TMP_PATH + 'temp_pepstats.pepstats')
    pepstats_dic_aas = pepstats(DATASET, TMP_PATH + 'temp_pepstats.pepstats')

    return pepstats_dic_aas
# -----------------------------------------------------------------------------------------------------------
def write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES):
    """ Function: write_FASTA_short_ids()

        Purpose:  Given a list of identifiers and the corresponding list 
                  of sequence, write these to a FASTA file using short
                  identifiers such as protein1, protein2, .... This is 
                  done because some programs like pepstats do not like 
                  long identifier names as input.
              
        Input:    Path to desired FASTA format output file, list of 
                  identifiers and list of corresponding sequences.
    
        Return:   List of short identifiers.
    """

    with open(f_output, 'w') as f:
        SHORT_IDENTIFIERS = []
        # Change identifiers to protein1, protein2, ...
        # and write to temporary file
        SET = zip(ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES)
    
        for index,  (identifier, sequence) in enumerate(SET):
            short_id = '>protein' + str(index)
            SHORT_IDENTIFIERS.append(short_id)
            f.writelines(short_id + '\n')
            f.writelines(sequence + '\n')

    return SHORT_IDENTIFIERS
# -----------------------------------------------------------------------------------------------------------
def pepstats(DATASET, FILE):
   
    pepstats_dic = {}    

    pepstats_file = FILE

    f = open(pepstats_file, 'r')
    content = f.readlines()
    f.close()
    
    for start, line in enumerate(content):
        if 'PEPSTATS of ' in line and ' from 1 to ' in line:		
            TARGET_ID = line.split('PEPSTATS of ')[1]
            TARGET_ID = TARGET_ID.split(' from')[0].strip()
            length = line.split(' to ')[1].strip()
            amino_acid_frequencies, dayhoff = [], []
            aminoacids_line = content[start + 10:start + 36]

            for item in aminoacids_line:	
                # Residue		Number		Mole%		DayhoffStat	    
                if item.split('\t')[0] != 'B = Asx' and item.split('\t')[0] != 'J = ---' and item.split('\t')[0] != 'O = ---' and item.split('\t')[0] != 'U = ---' and item.split('\t')[0] != 'Z = Glx' and item.split('\t')[0] != 'X = Xaa':
                    amino_acid_frequencies.append(float(item.split('\t')[4]))
                    dayhoff.append(float(item.split('\t')[6]))			
            # Extract molecular weight
            molecular_weight_line = content[start + 2:start + 3]
            molecular_weight = float(re.findall("\d+.\d+", str(molecular_weight_line))[0])
            # Extract charge
            charge_line = content[start + 3:start + 4]
            charge = float(re.findall("[-+]?\d+.\d+", str(charge_line))[1])

            # Extract isoelectric point
            isoelectric_line = content[start + 4:start + 5]
            isoelectric = float(re.findall("\d+.\d+", str(isoelectric_line))[0])

		# Watch out, in the pepstats software, if isoelectric point == None, an 
                # extra line will be introduced
            start_aas = content[start:].index('Property\tResidues\t\tNumber\t\tMole%\n')
            perline = content[start + start_aas + 1:start + start_aas + 10]

            tiny = float(re.findall("\d+.\d+", str(perline[0]))[-1])
            small = float(re.findall("\d+.\d+", str(perline[1]))[-1])
            aliphatic = float(re.findall("\d+.\d+", str(perline[2]))[-1])
            aromatic = float(re.findall("\d+.\d+", str(perline[3]))[-1])
            non_polar = float(re.findall("\d+.\d+", str(perline[4]))[-1])
            polar = float(re.findall("\d+.\d+", str(perline[5]))[-1])
            charged = float(re.findall("\d+.\d+", str(perline[6]))[-1])
            basic = float(re.findall("\d+.\d+", str(perline[7]))[-1])
            acidic = float(re.findall("\d+.\d+", str(perline[8]))[-1])
                    
            # Store values in dictionary
            pepstats_dic[TARGET_ID] = molecular_weight, charge, isoelectric, tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic, amino_acid_frequencies, length, dayhoff
                
    return pepstats_dic

def pepstats_old(DATASET, FILE):
   
    pepstats_dic = {}    

    pepstats_file = FILE

    f = open(pepstats_file, 'r')
    content = f.readlines()
    f.close()
    
    for (TARGET_ID, sequence) in DATASET:
        TARGET_ID = TARGET_ID.split(';')[0]
        TARGET_ID = TARGET_ID.replace('>','')
        TARGET_ID = TARGET_ID.strip()
        for start, line in enumerate(content):
            if 'PEPSTATS of ' + TARGET_ID + ' from 1 to ' in line:		
                length = float(len(sequence.strip()))
                amino_acid_frequencies, dayhoff = [], []
                aminoacids_line = content[start + 10:start + 36]

                for item in aminoacids_line:	
                    # Residue		Number		Mole%		DayhoffStat	    
                    if item.split('\t')[0] != 'B = Asx' and item.split('\t')[0] != 'J = ---' and item.split('\t')[0] != 'O = ---' and item.split('\t')[0] != 'U = ---' and item.split('\t')[0] != 'Z = Glx' and item.split('\t')[0] != 'X = Xaa':
                        amino_acid_frequencies.append(float(item.split('\t')[4]))
                        dayhoff.append(float(item.split('\t')[6]))			
                # Extract molecular weight
                molecular_weight_line = content[start + 2:start + 3]
                molecular_weight = float(re.findall("\d+.\d+", str(molecular_weight_line))[0])
                # Extract charge
                charge_line = content[start + 3:start + 4]
                charge = float(re.findall("[-+]?\d+.\d+", str(charge_line))[1])

                # Extract isoelectric point
                isoelectric_line = content[start + 4:start + 5]
                isoelectric = float(re.findall("\d+.\d+", str(isoelectric_line))[0])

		# Watch out, in the pepstats software, if isoelectric point == None, an 
                # extra line will be introduced
                start_aas = content[start:].index('Property\tResidues\t\tNumber\t\tMole%\n')
                perline = content[start + start_aas + 1:start + start_aas + 10]

                tiny = float(re.findall("\d+.\d+", str(perline[0]))[-1])
                small = float(re.findall("\d+.\d+", str(perline[1]))[-1])
                aliphatic = float(re.findall("\d+.\d+", str(perline[2]))[-1])
                aromatic = float(re.findall("\d+.\d+", str(perline[3]))[-1])
                non_polar = float(re.findall("\d+.\d+", str(perline[4]))[-1])
                polar = float(re.findall("\d+.\d+", str(perline[5]))[-1])
                charged = float(re.findall("\d+.\d+", str(perline[6]))[-1])
                basic = float(re.findall("\d+.\d+", str(perline[7]))[-1])
                acidic = float(re.findall("\d+.\d+", str(perline[8]))[-1])
                    
                # Store values in dictionary
                pepstats_dic[TARGET_ID] = molecular_weight, charge, isoelectric, tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic, amino_acid_frequencies, length, dayhoff
                
    return pepstats_dic
# -----------------------------------------------------------------------------------------------------------
def get_seqs_ids_fasta(FASTA_FILE):
    """ Function: get_seqs_ids_fasta()

        Purpose:  Given a FASTA format file, this function extracts
                  the list of identifiers and the list of sequences 
                  in the order in which they appear in the FASTA file.
              
        Input:    Path to FASTA format file.
    
        Return:   List of identifiers and list of sequences in the order 
                  in which they appear in the FASTA file.
    """ 
    identifiers = []
    sequences = []

    with open(FASTA_FILE) as f: 
        content = f.readlines()

        for position, line in enumerate(content):
            if '>' in line:
                identifiers.append(line)
                seq = []
                following_lines = content[position + 1:]
                for next_line in following_lines:
                    if '>'not in next_line:
                        seq.append(next_line.strip())
                    else:
                        break
                sequence = "".join(seq)
                sequence = sequence.replace('*', '')
                sequences.append(sequence)

    return identifiers, sequences
# -----------------------------------------------------------------------------------------------------------
def seq_homology_testset(identifier, sequence, TESTSET):

    f = open('seq_test.fasta', 'w')
    f.writelines(identifier)
    f.writelines(sequence + '\n')
    f.close()

    f = open('testset.fasta', 'w')
    for ident, seq in TESTSET:         
        f.writelines(ident)
        f.writelines(seq + '\n')
    f.close()

    #print "Scan with phmmer"
    os.popen('phmmer -E 0.00001 --tblout phmmer_seq_test.temp seq_test.fasta testset.fasta')
    f = open('phmmer_seq_test.temp', 'r')
    content = f.readlines()
    f.close()

    homologs = []
    HOMOLOG_IN_SAMPLE = False
    for line in content:   
        if not '#' in line:        
            homolog = line.split()[0]
            homolog = '>' + homolog + '\n'
            homologs.append(homolog)
            if homologs and len(homologs) > 1:
                HOMOLOG_IN_SAMPLE = True
    
    if not HOMOLOG_IN_SAMPLE and not (identifier, sequence.replace('*','')) in TESTSET:
        TESTSET.append((identifier, sequence.replace('*','')))

    return TESTSET
# -----------------------------------------------------------------------------------------------------------
def filterX(IDENTIFIERS, SEQUENCES):
    # Replace ambiguous amino acids because ProtParam can't deal with them later on 
    # B = Asx = Aspartic acid or Asparagine
    # Z = Glx = Glutamic acid or Glutamine
    # X = any amino acid
    replaceB = ['D', 'N']
    replaceZ = ['E', 'Q']
    replaceU = ['C']
    replaceX = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    IDENTIFIERS_FILTERED, SEQUENCES_FILTERED = [], []
    
    for identifier, sequence in zip(IDENTIFIERS, SEQUENCES):

        if 'B' in sequence or 'Z' in sequence or 'X' in sequence or 'U' in sequence:
            IDENTIFIERS_FILTERED.append(identifier)

            sequence_replaced = []
            for char in sequence:
                char = char.replace('B', random.choice(replaceB))
                char = char.replace('Z', random.choice(replaceZ))
                char = char.replace('X', random.choice(replaceX))
                char = char.replace('U', random.choice(replaceU))
                sequence_replaced += char
            sequence_replaced = "".join(sequence_replaced)

            SEQUENCES_FILTERED.append(sequence_replaced.replace('*',''))
        else:
            IDENTIFIERS_FILTERED.append(identifier)
            SEQUENCES_FILTERED.append(sequence.replace('*',''))

    return IDENTIFIERS_FILTERED, SEQUENCES_FILTERED 
# -----------------------------------------------------------------------------------------------------------
def apoplast_classifier(pepstats_dic_aas, TARGET_IDs, sequences, TMP_PATH, MODEL):

    predicted_nucleus, predictions = [], []
    weka = TMP_PATH + 'apoplast.arff'
    X = [[] for __ in xrange(len(TARGET_IDs))]

    protein_position = 0

    with open(weka, 'w') as f:

        for TARGET_ID, sequence in zip(TARGET_IDs, sequences):

            TARGET_ID = TARGET_ID.replace('>', '')
            TARGET_ID = TARGET_ID.strip()
            # Create a list of features for each protein
            
            molecular_weight, charge, isoelectric, tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic, amino_acid_frequencies, length, dayhoff = pepstats_dic_aas[TARGET_ID]

            prot = ProtParam.ProteinAnalysis(sequence.replace('*',''))            
            
            #X[protein_position] = [charge, isoelectric] + [tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic] + amino_acid_frequencies + [GRAVY(sequence)] + [prot.aromaticity(), prot.instability_index(), prot.molecular_weight()] + [sequence.count('C')] 
            X[protein_position] = [charge, isoelectric] + [tiny, small, aliphatic, aromatic, non_polar, polar, charged, basic, acidic] + amino_acid_frequencies + [GRAVY(sequence)] + [prot.aromaticity(), prot.instability_index()] + [sequence.count('C')] 

            protein_position += 1            

        f.writelines(ARFF_APOPLAST_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    ParamList = ['java', '-cp', '/home/spe12g/Work_Csiro/Software/weka-3-8-1/weka.jar', 'weka.classifiers.trees.RandomForest',  
             '-l', MODEL, '-T', weka, '-p', 'first-last']

    with open(TMP_PATH + 'APOPLAST_Predictions.txt', 'wb') as out:
        Process = subprocess.Popen(ParamList, shell=False, stdout=out)
        sts = Process.wait()

    return 
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def short_output(predictions):
    """ Function: short_output()

        Purpose:  Given the WEKA predictions for each protein, write  
                  string that contains the short output format.
              
        Input:    WEKA predictions for each protein.                  
    
        Return:   String that contains predictions for all proteins as tab-delimited table.
    """
    # Output predictions for all proteins as tab-delimited table
    short_output_string = '# Identifier \t Prediction \t Probability \n'
    for protein, pred, prob, sequence in predictions:    
        short_output_string += protein + '\t' + pred + '\t' + str(prob) + '\n'            

    return short_output_string
# -----------------------------------------------------------------------------------------------------------
def long_output(ORIGINAL_IDENTIFIERS, predicted_effectors):
    """ Function: long_output()

        Purpose:  Given the predicted effectors and identifiers for the test set,  
                  write string that contains the long output format.
              
        Input:    Predicted effectors and identifiers of test set.                  
    
        Return:   String that contains list of predicted effectors with posterior probabilites
                  and a short statistic on the percentage of predicted effectors in the test set.
    """
    # Output predicted effectors for long format
    long_output_string = '-----------------\n'
    long_output_string += 'Predicted apoplastic proteins:\n\n'
    for effector, prob, sequence in predicted_effectors:
        long_output_string += effector + '| Apoplastic probability:' + str(prob) + '\n'

    long_output_string += '-----------------\n\n'
    long_output_string += 'Number of proteins that were tested: ' + str(len(ORIGINAL_IDENTIFIERS)) + '\n' 
    long_output_string += 'Number of predicted apoplastic proteins: ' + str(len(predicted_effectors)) + '\n' 
    long_output_string += '\n' + '-----------------' + '\n' 
    long_output_string += str(round(100.0*len(predicted_effectors)/len(ORIGINAL_IDENTIFIERS), 1)) + ' percent are predicted to be apoplastic.'  
    long_output_string += '\n' + '-----------------' + '\n'

    return long_output_string
# -----------------------------------------------------------------------------------------------------------           

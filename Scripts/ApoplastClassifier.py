import os
import sys
import functions
import subprocess
import errno
import uuid
import shutil
import re
from operator import itemgetter
import numpy as np
import tempfile
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Main Program starts here
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    SCRIPT_PATH = sys.path[0]
    #MODEL = SCRIPT_PATH + '/ITERATIONS/RATIO4_69_MODEL.model'
    #MODEL = SCRIPT_PATH + '/ITERATIONS_APO/RATIO5_93_MODEL.model'   #Chosen by maximizing AUC
    MODEL = SCRIPT_PATH + '/ITERATIONS_APO_MWremoved/RATIO4_55_MODEL.model'   #Chosen by maximizing AUC
    # -----------------------------------------------------------------------------------------------------------
    FASTA_FILE = sys.argv[1]
    #output_file, short_format, apoplast_output = None, None, None#sys.argv[2]
    commandline = sys.argv[1:]
    # -----------------------------------------------------------------------------------------------------------
    if commandline:
        FASTA_FILE, short_format, output_file, apoplast_output = functions.scan_arguments(commandline)
	# If no FASTA file was provided with the -i option
        if not FASTA_FILE:
            print
            print 'Please specify a FASTA input file using the -i option!'
            functions.usage()
    else:
        functions.usage()
    # -----------------------------------------------------------------------------------------------------------
    # Temporary folder name identifier that will be used to store results
    FOLDER_IDENTIFIER = str(uuid.uuid4())
    RESULTS_PATH = tempfile.mkdtemp() + '/'
    print FASTA_FILE
    # -----------------------------------------------------------------------------------------------------------
    ORIGINAL_IDENTIFIERS, SEQUENCES = functions.get_seqs_ids_fasta(FASTA_FILE)
    SEQUENCES = [seq.upper() for seq in SEQUENCES]
    #SEQUENCES = [seq[20:] for seq in SEQUENCES]
    print 'Proteins to classify: ', len(ORIGINAL_IDENTIFIERS)
    PROTEINS_CLASSIFIED = len(ORIGINAL_IDENTIFIERS)
    # -----------------------------------------------------------------------------------------------------------
    # Replace ambiguous amino acids and filter short proteins with length < 30 aas
    ORIGINAL_IDENTIFIERS, SEQUENCES = functions.filterX(ORIGINAL_IDENTIFIERS, SEQUENCES)
    print 'After filtering for sequences with X and other unusual bases: ', len(ORIGINAL_IDENTIFIERS)
    # -----------------------------------------------------------------------------------------------------------
    # Write new FASTA file with short identifiers because pepstats can't handle long names
    f_output = RESULTS_PATH + 'short_ids.fasta'
    SHORT_IDENTIFIERS = functions.write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    # Call pepstats
    print 'Call pepstats...'
    ProcessExe = 'pepstats'
    ParamList = [ProcessExe, '-sequence', f_output, 
              '-outfile', RESULTS_PATH + FOLDER_IDENTIFIER + '.pepstats']
    try:
        Process = subprocess.Popen(ParamList, shell=False)
        sts = Process.wait()
        cstdout, cstderr = Process.communicate()

        if Process.returncode:
            raise Exception("Calling pepstats returned %s"%Process.returncode)
        if cstdout:
            pass
        elif cstderr:
            sys.exit()
    except:
        e = sys.exc_info()[1]
        print "Error calling pepstats: %s" % e
        sys.exit()
    print 'Done.'
    print
    # -----------------------------------------------------------------------------------------------------------
    # Parse pepstats file
    print 'Scan pepstats file'
    pepstats_dic_aas = functions.pepstats(zip(SHORT_IDENTIFIERS, SEQUENCES), RESULTS_PATH + FOLDER_IDENTIFIER + '.pepstats')
    print 'Done.'
    print
    # -----------------------------------------------------------------------------------------------------------   
    # -----------------------------------------------------------------------------------------------------------
    functions.apoplast_classifier(pepstats_dic_aas, SHORT_IDENTIFIERS, SEQUENCES, RESULTS_PATH, MODEL)
    file_output = RESULTS_PATH + 'APOPLAST_Predictions.txt'
    predicted_apoplast, not_predicted_apoplast, predictions = functions.parse_weka_output(file_output, ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------
    # If user wants the stdout output directed to a specified file
    if output_file:

        with open(output_file, 'wb') as out:
            # Short format: output predictions for all proteins as tab-delimited table
            if short_format:
                out.writelines(functions.short_output(predictions))
            # If the user wants to see the long format, output additional information and stats
            else:
                out.writelines(functions.short_output(predictions))
                out.writelines(functions.long_output(ORIGINAL_IDENTIFIERS, predicted_apoplast))                
        print 'Results were saved to output file:', output_file  

    else:
        # Short format: output predictions for all proteins as tab-delimited table to stdout
        if short_format:
            print functions.short_output(predictions)
        # If the user wants to see the long format, output additional information and stats
        else:
            print functions.short_output(predictions)
            print functions.long_output(ORIGINAL_IDENTIFIERS, predicted_apoplast)
    # -----------------------------------------------------------------------------------------------------------
    # If the user additionally wants to save the predicted apoplastic proteins in a provided FASTA file
    if apoplast_output:
        with open(apoplast_output + '_apoplastic.fasta', 'w') as f_output:
            for apoplast, prob, sequence in predicted_apoplast:
                f_output.writelines('>' + apoplast + ' | Apoplast probability: ' + str(prob) + '\n')
                f_output.writelines(sequence + '\n')  
   
        with open(apoplast_output + '_non_apoplastic.fasta', 'w') as f_output:
            for apoplast, prob, sequence in not_predicted_apoplast:
                f_output.writelines('>' + apoplast + ' | Non-apoplast probability: ' + str(prob) + '\n')
                f_output.writelines(sequence + '\n')  
    # -----------------------------------------------------------------------------------------------------------
    # Clean up and delete temporary folder that was created
    shutil.rmtree(RESULTS_PATH)
    # -----------------------------------------------------------------------------------------------------------
    average_prob = sum([item[1] for item in predicted_apoplast])/len(predicted_apoplast)
    print 'Average probability:', round(average_prob,2)
    # -----------------------------------------------------------------------------------------------------------


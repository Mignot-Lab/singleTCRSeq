"""
scripts to process TCR fasta sequencing file to derive junctinonal indices
@author: Aditya Ambati ambati@stanford.edu, Mignot Lab, Stanford University
"""

import argparse
import csv
import logging
import os
import re
import subprocess
from collections import defaultdict
from datetime import datetime
from subprocess import PIPE


#colNeeded = ('Subject', 'Peptide', 'Dx', 'CDR3a', 'TCRa sequence', 'CDR3b', 'TCRb sequence', 'alt CDR3a', 'alt TCRa sequence')
#colNeeded = ('Subject', 'Peptide', 'Dx', 'CDR3a', 'TCRa.sequence', 'CDR3b', 'TCRb.sequence', 'alt.CDR3a', 'alt.TCRa.sequence')
def readData(filein, getCol=False, colNeeded = ('Subject', 'Peptide', 'Dx', 'CDR3a', 'TCRa.sequence', 'CDR3b', 'TCRb.sequence', 'alt.CDR3a', 'alt.TCRa.sequence')):
    '''reader function of single cell summary file'''
    data=[]
    headerDict =defaultdict(int)
    header = ''
    with open(filein) as sseq_in:
        #logging.info('READING DATA FROM FILE  {} '.format(filein))
        for n, line in enumerate(sseq_in):
            dataParse=line.strip('\n').split(',')
            if n > 0:
                    dataList = ','.join([dataParse[n] for n in colId]) +'$'+line.strip('\n')
                    data.append(dataList)
            else:
                headerDict={item.strip():n for n, item in enumerate(dataParse)}
                header += line.strip('\n')
                data.append(header)
                colId = [headerDict.get(item) for item in colNeeded]
                if getCol is True:
                    return colId
    return data

    


def igBlast(nucFasta, headFasta):
    '''Main blast engine repeatedly is called for every cdr3 fasta to make blast queries, expects a nucleotide file and a header string and returns parsed stdout from blast'''
    getDir = os.getcwd()
    if os.path.exists(getDir+'/temp/'+'fastaQueries'):
        fileFasta = getDir+'/temp/fastaQueries/'+headFasta+'.fasta'
    else:
        os.mkdir(getDir+'/temp/'+'fastaQueries')
        fileFasta = getDir+'/temp/fastaQueries/'+headFasta+'.fasta'

    with open(fileFasta, 'w') as outFasta:
        outFasta.write('>'+headFasta+'\n')
        outFasta.write(nucFasta+'\n')
    igCall = 'blast/bin/./igblastn -germline_db_V blast/database/Hu_V_processed.fasta -germline_db_J blast/database/Hu_J_processed.fasta -germline_db_D blast/database/Hu_D_processed.fasta -organism human -ig_seqtype TCR -num_alignments_V 1 -num_alignments_J 1 -num_alignments_D 1 -auxiliary_data blast/optional_file/human_gl_j24.aux -show_translation -outfmt 19 -query '+fileFasta
    print(igCall)
    igCallsub = subprocess.Popen(igCall, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = igCallsub.communicate()
    if not stderr:
        os.remove(fileFasta)
        return stdout
    else:
        logging.info('BLAST FAILED {}'.format(fileFasta))
        raise ValueError('Blast Not Possible {}'.format(fileFasta))


def parseBlast(stdOut, chain):
    '''expects the blast stdout and the chain of the cdr3sequence and returns a string with junctional indices for alpha and beta chains
    '''
    parseCall = stdOut.split('\t')
    seqFasta = parseCall[1]
    #vStart = int(parseCall[60])
    vEnd = int(parseCall[61])
    jStart = int(parseCall[68])
    #jEnd = int(parseCall[69])
    cdr3St = int(parseCall[82])-3
    cdr3En = parseCall[83]
    cdr3 = parseCall[43]
    if chain == 'beta':
        if parseCall[64]: # check if a d match is present
            dStart = int(parseCall[64])
            dEnd = int(parseCall[65])
            vdJunction = seqFasta[vEnd:dStart-1]
            djJunction = seqFasta[dEnd:jStart-1] 
            # make the string in the format <CDR3NUCIndex,vdjunctionNuc:djJunctionNuc,VendIndex,DstartIndex,DendIndex,jStart>
            makeStr = ' '.join([str(cdr3St)+':'+str(cdr3En), vdJunction+':'+djJunction, str(vEnd), str(dStart), str(dEnd), str(jStart)])
        else: # assume that no d matches and make a vj junction
            logging.info('D matches not present in  {}'.format(cdr3))
            vjJunction=seqFasta[vEnd:jStart-1]
            makeStr = ' '.join([str(cdr3St)+':'+str(cdr3En), 'NA:'+vjJunction, str(vEnd), 'NA', 'NA', str(jStart)])
    else:# alphas have only the vj so we are golden
        # make the string in the format <CDR3NUCIndex,vjJunctionNuc,VendIndex,jStart>
        vjJunction = seqFasta[vEnd:jStart-1]
        makeStr = ' '.join([str(cdr3St)+':'+str(cdr3En), vjJunction, str(vEnd), str(jStart)])
    return makeStr

def processData(fileIn):
    ''' reads the datalist item and process the cdr3 to derive junctions which are written to temp files and 
    also returned as dicts that can be input into parse junction funct
    '''
    data = readData(filein=fileIn) # length of 2 data[0] contains the calls, data[1] contains the header
    #colNeeded = ('Subject', 'Peptide', 'Dx', 'CDR3a', 'TCRa sequence', 'CDR3b', 'TCRb sequence', 'alt CDR3a', 'alt TCRa sequence')
    alphaList = set()
    betaList = set()
    alphaOut = open('temp/AlphaJunctions.csv', 'w')
    betaOut = open('temp/BetaJunctions.csv', 'w')
    alphasProcessed, betasProcessed = 0,0
    alphaJunc = {}
    betaJunc = {}
    # For every instance
    for n,line in enumerate(data):
        if n > 0:
            lineParse = line.split('$')
            parse_line=lineParse[0].strip('\n').split(',')
            if parse_line[0] and parse_line[1] and parse_line[2]:# only if sample id is present must we process
                make_key = ','.join([parse_line[0].strip(), parse_line[1].strip(), parse_line[2].strip()])
            cdr3a = parse_line[3]
            cdr3b = parse_line[5]
            cdr3a_alt= parse_line[7]
            alphaNuc = parse_line[4]
            betaNuc = parse_line[6]
            alphaAltNuc = parse_line[8]
            if betaNuc:
                if betaNuc not in betaList:
                    betaList.add(betaNuc)
                    blastCall=igBlast(nucFasta = betaNuc, headFasta=make_key+'_'+cdr3b)
                    if blastCall:
                        betasProcessed += 1
                        parseBlastcall = blastCall.decode().strip().split('\n')[-1]
                        outStr=parseBlast(stdOut=parseBlastcall, chain='beta')
                        betaOut.write(betaNuc+','+outStr+'\n')
                        betaJunc[betaNuc]=outStr
            if alphaNuc:
                if alphaNuc not in alphaList:
                    alphaList.add(alphaNuc)
                    blastCall=igBlast(nucFasta = alphaNuc, headFasta=make_key+'_'+cdr3a)
                    if blastCall:
                        alphasProcessed += 1
                        parseBlastcall = blastCall.decode().strip().split('\n')[-1]
                        outStr=parseBlast(stdOut=parseBlastcall, chain='alpha')
                        alphaOut.write(alphaNuc+','+outStr+'\n')
                        alphaJunc[alphaNuc]=outStr
            if alphaAltNuc:
                if alphaAltNuc not in alphaList:
                    alphaList.add(alphaAltNuc)
                    blastCall=igBlast(nucFasta = alphaAltNuc, headFasta=make_key+'_'+cdr3a_alt)
                    if blastCall:
                        alphasProcessed += 1
                        parseBlastcall = blastCall.decode().strip().split('\n')[-1]
                        outStr=parseBlast(stdOut=parseBlastcall, chain='alpha')
                        alphaOut.write(alphaAltNuc+','+outStr+'\n')
                        alphaJunc[alphaAltNuc]=outStr
    logging.info('PROCESSED {} BETAS AND {} ALPHA BLASTS'.format(betasProcessed, alphasProcessed))
    logging.info('JUNCTION FILES WRITTEN TO {} & {}'.format('AlphaJunctions.csv','BetaJunctions.csv'))
    alphaOut.close()
    betaOut.close()

def main():
    parser = argparse.ArgumentParser(description='A script to derive junctional indices from TCR fasta sequences')
    parser.add_argument('-TCRFile', help='Summary csv file of single cell TCR calls from HIMC stanford', required=True)
    args=parser.parse_args()
    fileIn = args.TCRFile
    now = datetime.now()
    logFile='SingleTCR_'+now.strftime("%m_%d_%Y_%H_%M_%S")+'.log'
    logging.basicConfig(filename=logFile,level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logging.info('INPUT FILE {}'.format(fileIn))
    processData(fileIn=fileIn)

if __name__ == "__main__":main()

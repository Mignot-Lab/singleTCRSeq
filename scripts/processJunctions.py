"""
scripts to process TCR fasta sequencing file to annotate junctinonal indices
@author: Aditya Ambati ambati@stanford.edu, Mignot Lab, Stanford University
"""
import argparse
import logging
from collections import defaultdict
from datetime import datetime
import processSingleSeq


def parseJuncfiles(fileIn):
    nucDict = defaultdict(str)
    with open(fileIn) as jFile:
        for line in jFile:
            lineParse = line.strip().split(',')
            nucKey = lineParse[0]
            valueDict = lineParse[1]#','.join(lineParse[1:])
            nucDict[nucKey] = valueDict
    return nucDict


def annotateCalls(fileIn, alphaJuncDict, betaJuncDict):
    '''This will take the single cell summary file and annotate every cdr3 in the string
    <CDR3NUCIndex,vdjunctionNuc:djJunctionNuc,VendIndex,DstartIndex,DendIndex,jStart> for beta
    while alpha is annotated <CDR3NUCIndex,vjJunctionNuc,VendIndex,jStart>
    '''
    ## instantiate an empty file to annotate
    outFile = open(fileIn.replace('.csv', 'JunctionAnnotated.csv'), 'w')
    ## read the data file
    data=processSingleSeq.readData(fileIn)
    #colVec = processSingleSeq.readData(fileIn, getCol=True)
    outFile.write(data[0]+',BetaJunc,AlphaJunc,AltAlphaJunc'+'\n')
    for n,line in enumerate(data):
        if n > 0:
            lineParse = line.split('$')
            outFile.write(lineParse[1]+',')
            parse_line=lineParse[0].strip('\n').split(',')
            alphaNuc = parse_line[4]
            betaNuc = parse_line[6]
            alphaAltNuc = parse_line[8]
            if betaNuc:
                if betaNuc in betaJuncDict:
                    getValBeta = betaJuncDict.get(betaNuc)
                    outFile.write(getValBeta+',')
            else:
                outFile.write(',')
            if alphaNuc:
                if alphaNuc in alphaJuncDict:
                    getValAlpha = alphaJuncDict.get(alphaNuc)
                    outFile.write(getValAlpha+',')
            else:
                outFile.write(',')
            if alphaAltNuc:
                if alphaAltNuc in alphaJuncDict:
                    getValAlphaAlt = alphaJuncDict.get(alphaAltNuc)
                    outFile.write(getValAlphaAlt+'\n')
            else:
                outFile.write(','+'\n')

# def main():
#     parser = argparse.ArgumentParser(description='A script to derive junctional indices from TCR fasta sequences')
#     parser.add_argument('-TCRFile', help='Summary csv file of single cell TCR calls from HIMC stanford', required=True)
#     args=parser.parse_args()
#     fileIn = args.TCRFile
#     now = datetime.now()
#     logFile='AnnotationJunction_'+now.strftime("%m_%d_%Y_%H_%M_%S")+'.log'
#     logging.basicConfig(filename=logFile,level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#     logging.info('INPUT FILE {}'.format(fileIn))
#     alphaJuncDict = parseJuncfiles(fileIn='temp/AlphaJunctions.csv')
#     betaJuncDict = parseJuncfiles(fileIn='temp/BetaJunctions.csv')
#     annotateCalls(fileIn=fileIn, alphaJuncDict=alphaJuncDict, betaJuncDict=betaJuncDict)

# if __name__ == "__main__":main()

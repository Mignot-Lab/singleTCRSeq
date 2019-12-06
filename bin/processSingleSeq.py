from collections import defaultdict
import subprocess
from subprocess import PIPE
import os
import re


filein = '/Users/ambati/Documents/SingleCellTCR/COMBINED_SSEQ_GUO_October12_2018_edit_fornorm_lowthreshold.csv'

def readData(filein):
    '''reader function of single cell summary file'''
    data=[]
    with open(filein) as sseq_in:
        print('READING DATA FROM FILE  {} '.format(filein))
        for line in sseq_in:
            data.append(line.strip('\n'))

    return data


def igBlast(nucFasta, headFasta):
    '''Function to make blast queries'''
    getDir = os.getcwd()
    if os.path.exists(getDir+'/'+'fastaQueries'):
        fileFasta = getDir+'/fastaQueries/'+headFasta+'.fasta'
    else:
        os.mkdir(getDir+'/'+'fastaQueries')
        fileFasta = getDir+'/fastaQueries/'+headFasta+'.fasta'

    with open(fileFasta, 'w') as outFasta:
        outFasta.write('>'+headFasta+'\n')
        outFasta.write(nucFasta+'\n')
    igCall = '/Users/ambati/ncbi-igblast-1.14.0/bin/./igblastn -germline_db_V /Users/ambati/ncbi-igblast-1.14.0/database/Hu_V_processed.fasta -germline_db_J /Users/ambati/ncbi-igblast-1.14.0/database/Hu_J_processed.fasta -germline_db_D /Users/ambati/ncbi-igblast-1.14.0/database/Hu_D_processed.fasta -organism human -ig_seqtype TCR -num_alignments_V 1 -num_alignments_J 1 -num_alignments_D 1 -auxiliary_data /Users/ambati/ncbi-igblast-1.14.0/optional_file/human_gl.aux -show_translation -outfmt 19 -query '+fileFasta
    print(igCall)
    igCallsub = subprocess.Popen(igCall, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = igCallsub.communicate()
    if not stderr:
        return stdout
    else:
        raise ValueError('Blast Not Possible {}'.format(fileFasta))


def parseBlast(stdOut, chain):
    parseCall = stdOut.split('\t')
    seqFasta = parseCall[1]
    vStart = int(parseCall[60])
    vEnd = int(parseCall[61])
    jStart = int(parseCall[68])
    jEnd = int(parseCall[69])
    cdr3 = parseCall[42]
    #cdr3St = int(parseCall[82])
    #cdr3En = int(parseCall[83])
    if chain == 'beta':
        try:
            dStart = int(parseCall[64])
            dEnd = int(parseCall[65])
            vdJunction = seqFasta[vEnd:dStart]
            djJunction = seqFasta[dEnd:jStart]
            makeStr = ','.join([str(vStart)+':'+str(jEnd)+'>'+vdJunction+':'+djJunction, str(vEnd), str(dStart), str(dEnd), str(jStart)])
            # makeStr = ','.join([str(cdr3St)+':'+str(cdr3En)+'>'+vdJunction+':'+djJunction, str(vEnd), str(dStart), str(dEnd), str(jStart)])
            return makeStr
        except:
            ValueError
            vjJunction=seqFasta[vEnd:jStart-1]
            makeStr = ','.join([str(vStart)+':'+str(jEnd)+'>'+vjJunction+':DJonly', str(vEnd), 'NA','NA', str(jStart-1)])
            # makeStr = ','.join([str(cdr3St)+':'+str(cdr3En)+'>'+vjJunction+':DJonly', str(vEnd), 'NA','NA', str(jStart-1)])
            return makeStr
    else:
        vjJunction = seqFasta[vEnd:jStart-1]
        makeStr = ','.join([str(vStart)+':'+str(jEnd)+'>'+vjJunction, str(vEnd), str(jStart-1)])
        # makeStr = ','.join([str(cdr3St)+':'+str(cdr3En)+'>'+vjJunction, str(vEnd), str(jStart-1)])
        return makeStr


dataList = readData(filein=filein)
alphaList = set()
betaList = set()
alphaOut = open('/Users/ambati/Documents/SingleCellTCR/AlphaJunctions.csv', 'w')
betaOut = open('/Users/ambati/Documents/SingleCellTCR/BetaJunctions.csv', 'w')
# For every instance
for n, line in enumerate(dataList):
    if n > 0:
        if ";" in line:
            parse_line=line.strip('\n').split(';')
        else:
            parse_line=line.strip('\n').split(',')
        if parse_line[1] and parse_line[2] and parse_line[3] and parse_line[0]:
            make_key = ','.join([parse_line[1].strip(), parse_line[2].strip(), parse_line[3].strip()])
        VB=parse_line[6]
        VA=parse_line[12]
        JA = parse_line[14].split(' ')[0]
        VAalt = parse_line[18]
        JAalt = parse_line[19].split(' ')[0]
        cdr3a = parse_line[15]
        cdr3b = parse_line[9]
        cdr3a_alt= parse_line[20]
        alphaNuc = parse_line[24]
        betaNuc = parse_line[23]
        alphaAltNuc = parse_line[25]
        if alphaNuc:
            if alphaNuc not in alphaList:
                alphaList.add(alphaNuc)
                blastCall=igBlast(nucFasta = alphaNuc, headFasta=make_key+'_'+cdr3a)
                if blastCall:
                    parseBlastcall = blastCall.decode().strip().split('\n')[-1]
                    outStr=parseBlast(stdOut=parseBlastcall, chain='alpha')
                    alphaOut.write(alphaNuc+','+outStr+'\n')
        if betaNuc:
            if betaNuc not in betaList:
                betaList.add(betaNuc)
                blastCall=igBlast(nucFasta = betaNuc, headFasta=make_key+'_'+cdr3b)
                if blastCall:
                    parseBlastcall = blastCall.decode().strip().split('\n')[-1]
                    outStr=parseBlast(stdOut=parseBlastcall, chain='beta')
                    betaOut.write(betaNuc+','+outStr+'\n')
        if alphaAltNuc:
            if alphaAltNuc not in alphaList:
                alphaList.add(alphaAltNuc)
                blastCall=igBlast(nucFasta = alphaAltNuc, headFasta=make_key+'_'+cdr3a_alt)
                if blastCall:
                    parseBlastcall = blastCall.decode().strip().split('\n')[-1]
                    outStr=parseBlast(stdOut=parseBlastcall, chain='alpha')
                    alphaOut.write(alphaAltNuc+','+outStr+'\n')
alphaOut.close()            
betaOut.close() 


def parseJuncfiles(fileIn):
    nucDict = defaultdict(str)
    with open(fileIn) as jFile:
        for line in jFile:
            lineParse = line.strip().split(',')
            nucKey = lineParse[0]
            valueDict = ','.join(lineParse[1:])
            nucDict[nucKey] = valueDict
    return nucDict

alphaJuncDict = parseJuncfiles(fileIn = '/Users/ambati/Documents/SingleCellTCR/AlphaJunctions.csv')
betaJuncDict = parseJuncfiles(fileIn = '/Users/ambati/Documents/SingleCellTCR/BetaJunctions.csv')

platePos=list(map(lambda x: x.split(',')[4].replace('-', 'plate').replace('PL', 'plate'), dataList))
dataHeaders = dataList[0]
dataList.pop(0)
platePos.pop(0)
matchPlate=re.compile("^plate[0-9]+")
plateMaster = [re.findall(matchPlate, i)[0] for i in platePos]
if len(plateMaster) == len(dataList):
    print('Row Matched')

plateSpecificDict = defaultdict(list)
for plate, row in zip(plateMaster, dataList):
    plateSpecificDict[plate].append(row)

plateCalls = plateSpecificDict.get('plate14')
test =defaultdict(lambda:defaultdict(lambda:defaultdict(int)))
test =defaultdict(lambda:defaultdict(int))

n=0
platePos = set()
for line in plateCalls:
    n += 1
    if ";" in line:
        parse_line=line.strip('\n').split(';')
    else:
        parse_line=line.strip('\n').split(',')
    if parse_line[1] and parse_line[2] and parse_line[3] and parse_line[0]:
        make_key = ','.join([parse_line[1].strip(), parse_line[2].strip(), parse_line[3].strip()])
    platePos.add(parse_line[4])
    VB=parse_line[6]
    VA=parse_line[12]
    JA = parse_line[14].split(' ')[0]
    VAalt = parse_line[18]
    JAalt = parse_line[19].split(' ')[0]
    cdr3a = parse_line[15]
    cdr3b = parse_line[9]
    cdr3a_alt= parse_line[20]
    alphaNuc = parse_line[24]
    betaNuc = parse_line[23]
    alphaAltNuc = parse_line[25]
    if alphaNuc in alphaJuncDict:
        print(alphaJuncDict.get(alphaNuc), n)
        getJunc = alphaJuncDict.get(alphaNuc)
        #test[cdr3a].append(make_key+'>'+alphaNuc)
        test[cdr3a][alphaNuc + '<>' +getJunc] += 1
        #test[alphaNuc].append(make_key+'>'+cdr3a)


def plateSpecfic(plateCalls):



# test=igBlast(nucFasta = alphaNuc, headFasta='testalpha')

# test=igBlast(nucFasta = betaNuc, headFasta='testbeta')
# test2=test.decode().split('\n')
test2= blastCall.decode().split('\n')
n = 0
for i, j in zip(test2[0].split('\t'), test2[1].split('\t')):
    n += 1
    print(i, j, 'index>>', n-1)



# parseBlast(stdOut=test2[1], chain='beta')
# parseBlast(stdOut=test2[1], chain='alpha')

# parseCall = test2[1].split('\t')
# seqFasta = parseCall[1]
# vStart = int(parseCall[60])
# vEnd = int(parseCall[61])
# jStart = int(parseCall[68])
# jEnd = int(parseCall[69])
# if chain == 'beta':
#     dStart = int(parseCall[64])
#     dEnd = int(parseCall[65])
#     vdJunction = seqFasta[vEnd-1:dStart-1]
#     djJunction = seqFasta[dEnd-1:jStart-1]
#     print(vdJunction, djJunction)
# else:
#     vjJunction = seqFasta[vEnd-1:jStart-1]
#     print(vjJunction)
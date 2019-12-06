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
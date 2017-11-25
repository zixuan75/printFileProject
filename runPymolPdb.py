import os, glob, sys
line = ''
total = len(sys.argv)
if total < 3 :
    print('Must have at least than two arguments')
    sys.exit()
lineOneByOne=sys.argv[1]
outFile=sys.argv[2]
dnaFile=sys.argv[3]
if dnaFile:
    os.system("pymol extractedRunChains/"+lineOneByOne+" "+dnaFile+" > "+outFile)
else:
    os.system("pymol extractedRunChains/"+lineOneByOne+" > "+outFile)
lineByComp=''
core=open("sin2Core.txt", "w")
compRow=''
buriedAtom2ExposedAtom=''
exposedAtom2BuriedAtom=''

buriedAtom2SurfaceAtom=''
surfaceAtom2BuriedAtom=''
lineOfData=0
if total == 5:
    compRow=sys.argv[4]
    lineOfData=1
elif total == 6:
    buriedAtom2ExposedAtom=sys.argv[5]
    lineOfData=2
elif total == 7:
    exposedAtom2BuriedAtom=sys.argv[6]
    lineOfData=3
elif total == 8:
    buriedAtom2SurfaceAtom=sys.argv[7]
    lineOfData=4
elif total == 9:
    surfaceAtom2BuriedAtom=sys.argv[8]
    lineOfData=5
else:
    lineOfData=0
noOfCoreLine=1
pdbid=''
coreLine=''
coreLine2=''
lineOfAid=''
line2Buried=''
lineToLine=0
factor=1
idOfFactor=1
pdbidForFH=0
lineWithPdbid=0
line31018=''
sizeOfCnt=0
cntOfLine=0
line2Core=1
idLine=1
line2CoreAtomList = list()
naturalProtein=False
bindingProtein=False
isDNA=False
def square(a, b):
    return a*b
lc2LineList={'6S':'S','1D':'D','3C':'C','7Z':'Z'}
lc2Line = lc2LineList['6S']
lc2Line2 = lc2LineList['1D']
lc2Line3 = lc2LineList['3C']
lc2Line4 = lc2LineList['7Z']
with open("pdbLineOfPdb.txt") as pFile:
    content=pFile.readlines()
    for txtLine in content:
        clean_string = ' '.join(txtLine.split())
        tokens = clean_string.split()
        coreLine = txtLine[6:13].strip()
        pdbid = txtLine[0:4].strip()
        yetSince = txtLine[13:19].strip()
        coreLine2 = txtLine[19:25].strip()
        line31018 = txtLine[19:25].strip()
        line2Buried = coreLine.strip()+' '+line31018
        line2Buried = line2Buried+str(1)
        lineOfAid = coreLine.strip()
        line2Core = line2Core+1
        if pdbid == "1a22":
            line = pdbid
            pdbidForFH = lineOneByOne
            print("crLine=" + lineOfAid)
            line2Buried2CoreAtom = list(line2Buried)
            line2CoreAtomList = line2Buried2CoreAtom
            print('line2Core='+str(line2Core)+'.0')
            noOfTokens = len(tokens)
            line2OnCore=line
            core.write(line2OnCore+' '+yetSince+' '+coreLine+' '+line31018+"0\n")
            with open(pdbidForFH) as sFile:
                content=sFile.readlines()
                for txtLineInCore in content:
                    clean_string2 = ' '.join(txtLineInCore.split())
                    tokens2 = clean_string2.split()
                    cntOfLine=cntOfLine+1
                    if txtLineInCore[21:22]=="N" or txtLineInCore[21:22]=="A":
                        bindingProtein=True
                    elif txtLineInCore[21:22]=="F" or txtLineInCore[21:22]=="J" or txtLineInCore[21:22]=="B" or txtLineInCore[21:22]=="C":
                        naturalProtein=True
                        if txtLineInCore[17:20]==" DA" or txtLineInCore[17:20]==" DG" or txtLineInCore[17:20]==" DC" or txtLineInCore[17:20]==" DT":
                            isDNA=True
            sFile.close()
            lineWithPdbid=1
#            print("lineOneByOne="+lineOneByOne)
            print("cntOfLine="+str(cntOfLine))
            core.write('------------------------\n')
            idLine = float(line31018 )
            idOneByOne = str(len(coreLine2) + 805.8)
            for i in range(1, noOfTokens ):
                noOfCoreLine = i
                lineToLine=lineToLine+1
                core.write('|'+str(noOfCoreLine)+".00 "+line31018+" "+str(lineToLine)+".00 "+coreLine+'|\n')
                if idOneByOne < lineOfAid: 
                    idOneByOne = lineOfAid
            core.write('------------------------\n')
            sizeOfCnt = cntOfLine+float(line31018)
            dnaCnt=float(yetSince)+cntOfLine
            print("dnaCnt="+str(dnaCnt))
            lineByComp = str(sizeOfCnt+2)
            print("lineByCom="+str(lineByComp))
            lineCnt = square(float(sizeOfCnt), 2)-(float(coreLine)+float(lineOfAid))
            core.write('lineCnt='+str(lineCnt)+'\n')
            print("idOneByOne="+str(idOneByOne))
            print("sizeOfCnt="+str(sizeOfCnt))
pFile.close()
core.close()
filedata=lc2LineList
print('secondary='+filedata['6S'])

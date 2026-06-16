import pandas as pd
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np

from HWE import LH_Distrib, midPvalue, nbHete, ExpectedHete, ExpecteHomAlt

from copy import deepcopy
from liftover import get_lifter


import os.path
from os import path
import math




##################################################################################################
##############------------------ PARSING VCFs AND MAKING CSVs functions --------------------------
################################################################################################## 

def VCFtoCSV_Exome(vcf_file,medList,quantityList,ethnList,sexList,p_valueList,csv_file,includeHeader,
                   listAGAIN,listBPHUNTER,listHGMD,listpLOF,FSTs_header,Intervals):    
    
    converter = get_lifter('hg38', 'hg19', one_based=True)
    
    listChr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
 'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
    

    header = ['chromo','position','ID','ref','alt','qual','filter','n_alt_alleles']
    allPop=[]

    for j in ethnList:
        for i in medList:
            for k in sexList:
                for l in quantityList:
                    header.append(fieldName(i,l,j,k))
                    allPop.append(fieldName(i,l,j,k))
                                        
    allPop[0]='AC';allPop[1]='AN';allPop[2]='AF';allPop[3]='nhomalt';  #instead of using AC_total, AN_total etc
    
    for i in p_valueList:
        header.append(i)  
              
        
#     header.append('DP')
    header.append('segdup')
    header.append('lcr')
#     header.append('decoy')
    
    header.append('consequence')
    header.append('impact')
    header.append('symbol')
    header.append('gene')
    header.append('feature_type')
    header.append('feature')
    header.append('biotype')
    header.append('EXON')
    header.append('INTRON')
    header.append('HGVSc')
    header.append('HGVSp')
    header.append('canonical')
    header.append('vepMisc')
    
    header.append('AB_in_hete_distrib')
    header.append('DP_distrib')
    header.append('DP_count_larger')
    
    header.append('InbreedingCoeff')  
    
    header.append('AB_between_40_and_60')
    header.append('AB_average')
    header.append('AB_standard_deviation')
    
    header.append('DP_average')
    header.append('DP_standard_deviation')
        
    header.append('AGAIN')
    header.append('BPHUNTER')
    header.append('HGMD')
    header.append('pLOF')
    
    header.extend(FSTs_header)
    
    header.append('hg19_chromo')
    header.append('hg19_position')
    
    header.append('UM')
    
    header.append('CADD_phred')
    
    
    data=[]
    
    vcfFile= open(vcf_file,'r')
    vcfLines= vcfFile.readlines()
    
    
    rr=0
    for l in vcfLines:
        

        #if rr%1000000==0:
        #    print(str(rr) + ' lines read at ' + str(datetime.now()))


        if l[0] != '#':
            dataLine=[]

            lineSplit=l.split()
            
            if lineSplit[0] in listChr:
            

                #---------------------------------- CHR, POS, ID, REF, ALT, QUAL, FITLER --------------------

                for i in range(7):
                    dataLine.append(lineSplit[i])


                info=lineSplit[7]

                #---------------------------------- n_alt_alleles -------------------------

                dataLine.append(valueOfField('n_alt_alleles',info))

    #             print(allPop)

    #             print(len(allPop))


                #---------------------------------- AC, AN, AF, nhomalt for each population --------------------

                for i in allPop:
                    dataLine.append(valueOfField(i,info))




                #---------------------------------- p-values for each population --------------------


                for i in range(len(p_valueList)):

                    n=str(dataLine[8+i*4+1])
                    if n.isdigit() == True:
                        n=int(n)/2
                    else:
                        n=0

                    nA=str(dataLine[8+i*4])
                    if nA.isdigit() == True:
                        nA=int(nA)
                    else:
                        nA=0

                    nAA=str(dataLine[8+i*4+3])
                    if nAA.isdigit() == True:
                        nAA=int(nAA)
                    else:
                        nAA=0

                    nAB=nA-2*nAA
                    if n==0 or nA==0:
                        dataLine.append('error')
                    else:
                        dataLine.append(midPvalue(LH_Distrib(n,nA),nAB))


                #----------------------------------segdup, lcr --------------------


#                 dataLine.append(valueOfField('DP',info))

                if info.find(';segdup;')==-1:
                    dataLine.append(0)
                else:
                    dataLine.append(1)

                if info.find(';lcr;')==-1:
                    dataLine.append(0)
                else:
                    dataLine.append(1)

#                 if info.find(';decoy;')==-1:
#                     dataLine.append(0)
#                 else:
#                     dataLine.append(1)



                 #---------------------------------- VP info --------------------


                VEPs=valueOfField('vep',info).split(',')      

                vepInfo=[]
                for n in range(11):
                    vepInfo.append('none')
                cano='none'
                vepMisc='none'
                
                HGVSp = 'none'

                nbMisc=0
                nbReg=0


                for nbV in range(len(VEPs)):

                    elements=VEPs[nbV].split('|')

                    if len(elements)<8:
                        if nbMisc==0:
                            vepMisc+=(VEPs[nbV])
                        else:
                            vepMisc+=(';'+VEPs[nbV])
                        nbMisc+=1

                    else:
                        for n in range(1,12):
                            if nbReg==0:
                                vepInfo[n-1]+=(elements[n])
                            else:
                                vepInfo[n-1]+=(';'+(elements[n]))


                        if elements[24]=='YES'and nbReg==0:
                            cano+='YES'
                        if elements[24]=='YES'and nbReg>0:
                            cano+=';'+'YES'
                        if elements[24]!='YES'and nbReg==0:
                            cano+='NO'
                        if elements[24]!='YES'and nbReg>0:
                            cano+=';'+'NO'
                            
                            
#                         if nbReg==0:
#                             HGVSp+=str(elements[11])
#                         if nbReg>0:
#                             HGVSp+=';'+str(elements[11])

                        nbReg+=1

                for n in range(len(vepInfo)):
                    if len(vepInfo[n])>4:
                        vepInfo[n]=vepInfo[n][4:len(vepInfo[n])] 
                        
#                 if len(HGVSp)>4:
#                     HGVSp=HGVSp[4:len(HGVSp)]

                if len(vepMisc)>4:
                    vepMisc=vepMisc[4:len(vepMisc)]
                if len(cano)>4:
                    cano=cano[4:len(cano)]

                for n in range(len(vepInfo)):
                    dataLine.append(vepInfo[n])
                    
                    
#                 dataLine.append(HGVSp)
                dataLine.append(cano)
                dataLine.append(vepMisc)


                #----------------------------------AB, DP, Inbreeding coeff --------------------

                dataLine.append(valueOfField('ab_hist_alt_bin_freq',info))
                dataLine.append(valueOfField('dp_hist_all_bin_freq',info))
                dataLine.append(valueOfField('dp_hist_all_n_larger',info))
                dataLine.append(valueOfField('InbreedingCoeff',info))




                #----------------------------------AB between 0.40 and 0.60, average and std --------------------
                AB_distr = valueOfField('ab_hist_alt_bin_freq',info).split('|')  
                for i in range(len(AB_distr)):
                    AB_distr[i]=int(AB_distr[i])

                if sum(AB_distr) > 0:
                    ABbetween40and60=sum(AB_distr[8:12])/sum(AB_distr)*100
                else:
                    ABbetween40and60 = 'no_hete'

                AB_List = []
                for i in range(len(AB_distr)):
                    for j in range(AB_distr[i]):
                        AB_List.append(0.025+i*0.05)

                AB_List_np = np.array(AB_List)

                dataLine.append(ABbetween40and60)
                dataLine.append(np.mean(AB_List_np))
                dataLine.append(np.std(AB_List_np))

                #----------------------------------DP average and std --------------------

                #to make runtime of std much faster, did an approximation of of std by
                #divising the nb in each bin by a 100

                DP_distr = valueOfField('dp_hist_all_bin_freq',info).split('|') 
                DP_distr.append(valueOfField('dp_hist_all_n_larger',info))
                for i in range(len(DP_distr)):
                    DP_distr[i]=int(DP_distr[i])


                totalNb = sum(DP_distr)
                totalDP = 0
                for i in range(len(DP_distr)):
                    totalDP += (2.5+i*5)*DP_distr[i]
                    
                DP_List = []
                for i in range(len(DP_distr)):
                    for j in range(round(DP_distr[i]/100)):
                        DP_List.append(2.5+i*5)

                DP_List_np = np.array(DP_List)


                if totalNb > 0:
                    dataLine.append(totalDP/totalNb)
                else:
                    dataLine.append('division_by_zero')

                dataLine.append(np.std(DP_List_np))


                 #---------------------------BPHUNTER and AGAIN--------------------------------

                if ([dataLine[0],dataLine[1],dataLine[3],dataLine[4]] in listAGAIN) == True:
                    dataLine.append(1)
                else:
                    dataLine.append(0)


                if ([dataLine[0],dataLine[1],dataLine[3],dataLine[4]] in listBPHUNTER) == True:
                    dataLine.append(1)
                else:
                    dataLine.append(0)

                #--------------------------- HGMD --------------------------------
                infoID = [dataLine[0],str(dataLine[1]),dataLine[3],dataLine[4]]


                listHGMD_chr = listHGMD[listChr.index(str(dataLine[0]))]

                HGMD_tags = []
                for i in listHGMD_chr:
                    if (infoID==i[0:4]):
                        HGMD_tags.append(i[4])

                if len(HGMD_tags) == 0:
                    dataLine.append('no')
                else:
                    if 'DM' in HGMD_tags:
                        dataLine.append('DM')
                    else:
                        dataLine.append(';'.join(HGMD_tags))


                #--------------------------- pLOF --------------------------------
                conseq=dataLine[50].split(';')

    #             print(conseq)

                conseqsList=[]
                for j in range(len(conseq)):
                    conseqsList.append(conseq[j].split('&'))

                indexToKeep=[]

                for j in range(len(conseqsList)):
                    if len(intersection(conseqsList[j],listpLOF)) > 0:
                        indexToKeep.append(j)

                if len(indexToKeep)>0:
                    dataLine.append('1')
                else:
                    dataLine.append('0')



                #--------------------------- FSTs --------------------------------

                for n in range(len(ethnList)):
                    for j in range(n+1,len(ethnList)):
                        AC_pop1_1 = int(dataLine[8+n*4])
                        AN_pop1 = int(dataLine[9+n*4])

                        AC_pop2_1 = int(dataLine[8+j*4])
                        AN_pop2 = int(dataLine[9+j*4])

                        if n == 0:
                            AC_pop1_1 = AC_pop1_1 - AC_pop2_1
                            AN_pop1 = AN_pop1 - AN_pop2

                        dataLine.append(FST(AC_pop1_1, AN_pop1, AC_pop2_1, AN_pop2))

                #---------------------------------- hg19 position --------------------

                position = int(dataLine[1])
                chromo = str(dataLine[0])

                position_info = converter[chromo][position]

                if len(position_info) == 1 and (len(position_info[0])==3):
                    position_chromo = position_info[0][0]
                    position_nb = position_info[0][1]

                    if (position_chromo in listChr):

                        dataLine.append(str(position_chromo))
                        dataLine.append(str(position_nb))

                    else:
                        dataLine.append('No_lift')
                        dataLine.append('No_lift')
                else:
                    dataLine.append('No_lift')
                    dataLine.append('No_lift')


    #             print(len(dataLine))


                #---------------------------------- checking UM --------------------


                if dataLine[105] != 'No_lift':

                    listUM_chr = Intervals[listChr.index(str(position_chromo))]

                    if binary_search_intervals(listUM_chr, position_nb) == True:
                        dataLine.append(1)
                    else:
                        dataLine.append(0)

                else:
                    dataLine.append(1)
                    
                    
                #----------------------------------CADD --------------------

                dataLine.append(valueOfField('cadd_phred',info))

                #-----------------------------------------------------------------

                data.append(dataLine)
            



    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        if includeHeader==1:
            writer.writerow(header)

        writer.writerows(data)
        
        
def csvEachPop_Exome(variantsCSV,medList,quantityList,ethnList,sexList,p_valueList,namefile):
    ByEthn_data = []
    for i in range(len(ethnList)):
        ByEthn_data.append([])
    #print("\n");print("STARTS OPENING VCF FILE ", datetime.now());print("\n")
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()
    #print("\n");print("STARTS LOOPING THROUGH VCF FILE ", datetime.now());print("\n")
    
    rr=0
    for l in csvLines:
        #if rr%500000==0:
        #    print(str(rr) + ' lines read at ' + str(datetime.now()))
        if rr>=0:
            info = l[:-1].split(",")
#             print(len(info))
            for i in range(len(ethnList)):
                data=[]

                data.extend(info[0:8])
                
                if rr == 0:
                    data.extend(['AC', 'AN', 'AF', 'nhomalt','pv_HWE'])
                else:
                    data.extend(info[8+4*i:8+4*(i+1)])
                    data.append(info[40+i])
                
                data.extend(info[8:len(info)])  
#                 print(len(data))
                                    
                ByEthn_data[i].append(data)
        rr=rr+1
        
    for i in range(len(ethnList)):
        with open(namefile+ethnList[i]+'.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(ByEthn_data[i])        
        
    
    
def add_MORE_Exome(variantsCSV,newCSV):           
    #print("\n");print("PROGRAM STARTS AT ", datetime.now());print("\n")
    data = []
    #print("\n");print("STARTS OPENING VCF FILE ", datetime.now());print("\n")
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()
    #print("\n");print("STARTS LOOPING THROUGH VCF FILE ", datetime.now());print("\n")    
    header=csvLines[0][:-1].split(",")   
    
    header.append('Observed_Hete')
    header.append('Expected_Hete')
    header.append('Excess_Hete')
    
    header.append('Observed_HomAlt')
    header.append('Expected_HomAlt')
    header.append('Excess_HomAlt')
    
    header.append('Observed_HomRef')
    header.append('Expected_HomRef')
    header.append('Excess_HomRef')
    
    header.append('CorrectInbreedingCoeff')    
    
    header.append('chi_HomMinor')
    header.append('chi_Hete')
    header.append('chi_HomMajor')
    
    
    data.append(header)
    
    ii=0
    for line in csvLines[1:len(csvLines)]:
        
#         print(ii)
#         if ii%100==0:
#             print(str(ii) + ' lines read at ' + str(datetime.now()))
            
        info = line[:-1].split(",")
        

        #-------------------------Expected vs Observed----------------------------
                
        nA=int(info[8])
        n=int(info[9])/2
        altHomo=int(info[11])
        nbHomRef=n-altHomo-nbHete(nA,altHomo)
        ExpectedHomRef = n - ExpectedHete(n,nA) - ExpecteHomAlt(n,nA)
        
                       
        info.append(nbHete(nA,altHomo))
        info.append(ExpectedHete(n,nA))
        if nbHete(nA,altHomo) > ExpectedHete(n,nA):
            info.append(1)
        else:
            info.append(0)                        
            
        info.append(altHomo)
        info.append(ExpecteHomAlt(n,nA))
        if altHomo > ExpecteHomAlt(n,nA):
            info.append(1)
        else:
            info.append(0)            
                   
        info.append(nbHomRef)
        info.append(ExpectedHomRef)
        if nbHomRef > ExpectedHomRef:
            info.append(1)
        else:
            info.append(0)
                       
           
        
        info.append(1-nbHete(nA,altHomo)/ExpectedHete(n,nA))
          
            
        
        if ExpectedHomRef > ExpecteHomAlt(n,nA):
            ExpectedHomMinor=ExpecteHomAlt(n,nA)
            ObservedHomMinor=altHomo
            
            ExpectedHomMajor=ExpectedHomRef
            ObservedHomMajor=nbHomRef            
        else:
            ExpectedHomMinor=ExpectedHomRef
            ObservedHomMinor=nbHomRef
            
            ExpectedHomMajor=ExpecteHomAlt(n,nA)
            ObservedHomMajor=altHomo
            
        if ExpectedHomMinor > 0:
            info.append(((ExpectedHomMinor-ObservedHomMinor)**2)/ExpectedHomMinor)
        else:
            info.append(float('inf'))
            
        info.append(((ExpectedHete(n,nA)-nbHete(nA,altHomo))**2)/ExpectedHete(n,nA))      
        info.append(((ExpectedHomMajor-ObservedHomMajor)**2)/ExpectedHomMajor)
        


            
        #----------------------------------------------------------------------------                        
        data.append(info)
        ii+=1
    
    #print("\n");print("STARTS WRITING CSV", datetime.now());print("\n")
    with open(newCSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
    #print("\n");print("PROGRAM ENDS AT ", datetime.now());print("\n")
    
    
          

def Filter_filter(variantsCSV,newCSV):
    
    
    #print("\n");print("PROGRAM STARTS AT ", datetime.now());print("\n")
    data = []
    #print("\n");print("STARTS OPENING VCF FILE ", datetime.now());print("\n")
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()
    #print("\n");print("STARTS LOOPING THROUGH VCF FILE ", datetime.now());print("\n")
    
    header=csvLines[0][:-1].split(",")    
    data.append(header)
    
    for line in csvLines[1:len(csvLines)]:
        #if i%1000000==0:
        #    print(str(i) + ' lines read at ' + str(datetime.now()))
            
        info = line[:-1].split(",")
        
        if (info[6]=='PASS') or (info[6]=='InbreedingCoeff'):
            data.append(info)
    
    #print("\n");print("STARTS WRITING CSV", datetime.now());print("\n")
    with open(newCSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
    #print("\n");print("PROGRAM ENDS AT ", datetime.now());print("\n")
    
    
    
def Filter_BI(variantsCSV,newCSV):
    
    
    #print("\n");print("PROGRAM STARTS AT ", datetime.now());print("\n")
    data = []
    #print("\n");print("STARTS OPENING VCF FILE ", datetime.now());print("\n")
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()
    #print("\n");print("STARTS LOOPING THROUGH VCF FILE ", datetime.now());print("\n")
    
    header=csvLines[0][:-1].split(",")    
    data.append(header)
    
    #print(len(vcfLines))
    
    i=1
    while i in range(len(csvLines)):
        
        #if i%1000000==0:
        #    print(str(i) + ' lines read at ' + str(datetime.now()))
            
        info = csvLines[i][:-1].split(",")
        
        currentIndex=i
        
        if int(info[7]) > 1 and i<len(csvLines)-1:
            
            listAF=[]
            listAF.append(float(info[10]))
            #print(i)
            infoNext=csvLines[i+1].split(",")
            while infoNext[1] == info[1] and i<len(csvLines)-1:
                listAF.append(float(infoNext[10]))
                i=i+1
                
                if i<len(csvLines)-1:
                    infoNext=csvLines[i+1].split(",")
            
            if len(listAF)==1:
                data.append(info)
            
            else:
                indMax = listAF.index(max(listAF))
                if (sum(listAF) - listAF[indMax])<0.001:
                    
                    dataLine=csvLines[currentIndex+indMax][:-1].split(",")
                    data.append(dataLine)
                    
        elif i == len(csvLines)-1:
            data.append(info)
            
        elif int(info[7]) == 1:
            data.append(info)
      
        i+=1
        
    #print("\n");print("STARTS WRITING CSV", datetime.now());print("\n")
    with open(newCSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
    #print("\n");print("PROGRAM ENDS AT ", datetime.now());print("\n")
    
    
def Filter_Auto_AC_CallRate(variantsCSV,threshold_AC,threshold_CR,ethnNumbAl,newCSV):
    #print("\n");print("PROGRAM STARTS AT ", datetime.now());print("\n")
    data = []
    #print("\n");print("STARTS OPENING VCF FILE ", datetime.now());print("\n")
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()
    #print("\n");print("STARTS LOOPING THROUGH VCF FILE ", datetime.now());print("\n")
    
    header=csvLines[0][:-1].split(",")    
    data.append(header)
    
    i=1
    for line in csvLines[1:len(csvLines)]:
        #if i%100000==0:
            #print(str(i) + ' lines read at ' + str(datetime.now()))
            
        info = line[:-1].split(",")
        
        ac=info[8]
        an=info[9]
        
        callRate = float(info[9])/ethnNumbAl
        
        if int(ac)>threshold_AC and int(ac)<int(an)-threshold_AC and callRate>threshold_CR:
            if (str(info[0])!='X') and (str(info[0])!='Y') and (str(info[0])!='chrX') and (str(info[0])!='chrY'):
                data.append(info)
                
        i+=1;
            
    #print("\n");print("STARTS WRITING CSV", datetime.now());print("\n")
    with open(newCSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
    #print("\n");print("PROGRAM ENDS AT ", datetime.now());print("\n")
    
    
def valueOfField(field,string):
    
    if (field =='AC') or (field =='ALGORITHMS'):

        beg=(string.find(field+'='))+len(field)+1

        if beg == len(field):
            return -1

        i=0
        while string[beg+i]!=';':
            i=i+1    
        end = beg+i

        return string[beg:end]
        
    elif field =='vep':
        beg=(string.find(';'+field+'='))+len(field)+2

        if beg == len(field)+1:
            return -1
        
        return string[beg:len(string)]
        
    else:
        beg=(string.find(';'+field+'='))+len(field)+2

        if beg == len(field)+1:
            return -1

        i=0
        while string[beg+i]!=';':
            i=i+1    
        end = beg+i

        return string[beg:end]
    
     
def fieldName(med, quant, ethn, sex):
    if (med=='') & (ethn=='') & (sex==''):
        field=quant
    elif (med=='') & (ethn==''):
        field=quant+'_'+sex
    elif (ethn=='') & (sex==''):
        field=med+'_'+quant
    elif (med=='') & (sex==''):
        field=quant+'_'+ethn
    elif (med==''):
        field=quant+'_'+ethn+'_'+sex
    elif (ethn==''):
        field=med+'_'+quant+'_'+sex
    elif (sex==''):
        field=med+'_'+quant+'_'+ethn
    else:
        field=med+'_'+quant+'_'+ethn+'_'+sex
        
    return field



def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def binary_search_intervals(intervals, position):
    
    left = 0
    right = len(intervals) - 1
    
    i = 0
    while left <= right:
        mid = left + (right - left) // 2
        start, end = intervals[mid]
        

        if start <= position <= end:
            return True
        elif position < start:
            right = mid - 1
        else:
            left = mid + 1
        i+=1
    
    return False



# Imported but unused function stubs (to prevent import errors)
def add_FSTs(variantsCSV,newCSV,FSTs_header,ethnList):
    """Function stub - imported but not used by notebook"""
    pass

def uniqueElements(l):
    """Function stub - imported but not used by notebook"""
    pass

def addUM(variantsCSV,newCSV,listUM,listChr):
    """Function stub - imported but not used by notebook"""
    pass

def Filter_UM(variantsCSV,newCSV):
    """Function stub - imported but not used by notebook"""
    pass

def IntervalsSum(intervals):
    """Function stub - imported but not used by notebook"""
    pass

def mergeOverlappingIntervals(arr):
    """Function stub - imported but not used by notebook"""
    pass

def Filter_UM_pop(variantsCSV,newCSV):
    """Function stub - imported but not used by notebook"""
    pass

def FST(AC_pop1_1, AN_pop1, AC_pop2_1, AN_pop2):
    
    if ((AC_pop1_1 == 0) and (AC_pop2_1 == 0)) or ((AN_pop1 == 0) or (AN_pop2 == 0)):
        return 'division_by_0'
    
    if (AC_pop1_1 == AN_pop1) and (AC_pop2_1 == AN_pop2):
        return 'AC=AN_both_pop'
    
    AF_pop1_1 = AC_pop1_1 / AN_pop1
    AF_pop1_2 = 1 - AF_pop1_1
    
    AF_pop2_1 = AC_pop2_1 / AN_pop2
    AF_pop2_2 = 1 - AF_pop2_1
    
    alpha_1 = 1 - AF_pop1_1**2 - AF_pop1_2**2
    alpha_2 = 1 - AF_pop2_1**2 - AF_pop2_2**2
    
    a = (1/2) * ( (AF_pop1_1-AF_pop2_1)**2 + (AF_pop1_2-AF_pop2_2)**2 )
    a = a - (AN_pop1+AN_pop2) * (AN_pop1*alpha_1+AN_pop2*alpha_2) / ( 4*AN_pop1*AN_pop2 * (AN_pop1+AN_pop2-1) )
    
    
    a_plus_b = (1/2) * ( (AF_pop1_1-AF_pop2_1)**2 + (AF_pop1_2-AF_pop2_2)**2 )
    a_plus_b = a_plus_b + (4*AN_pop1*AN_pop2 - AN_pop1 - AN_pop2) * (AN_pop1*alpha_1+AN_pop2*alpha_2) / ( 4*AN_pop1*AN_pop2 * (AN_pop1+AN_pop2-1) )
        
    return a / a_plus_b

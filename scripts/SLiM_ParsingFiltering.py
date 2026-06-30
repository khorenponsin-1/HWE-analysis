import pandas as pd
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib as mpl
import numpy as np

from HWE_slim import LH_Distrib, midPvalue, nbHete, ExpectedHete, ExpecteHomAlt, degreeOfDiseq, CE_Distrib, ExpecteHomAlt_CE, ExpectedHete_CE

import os.path
from os import path
import math


def find_all(string,seq):

    import re
    a_string = string
    indices_object = re.finditer(pattern=seq, string=a_string)
    indices = [index.start() for index in indices_object]
    return indices

def uniqueElements(l):
    output=[]
    for i in l:
        if i not in output:
            output.append(i)
    return output

def listIndOc(string, char):
    l=[]
    for i in range(len(string)):
        if string[i]==char:
            l.append(i)
    return l

def valueOfField(field,string):
    
    if field =='MID':

        beg=(string.find(field+'='))+len(field)+1

        if beg == len(field):
            return -1

        i=0
        while string[beg+i]!=';':
            i=i+1    
        end = beg+i

        return string[beg:end]
        
    elif field =='DP':
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


def populationPresent(med, ethn, sex, string):
    if (med=='') & (ethn=='') & (sex==''):
        field='AC'
    elif (med=='') & (ethn==''):
        field='AC'+'_'+sex
    elif (ethn=='') & (sex==''):
        field=med+'_'+'AC'
    elif (med=='') & (sex==''):
        field='AC'+'_'+ethn
    elif (med==''):
        field='AC'+'_'+ethn+'_'+sex
    elif (ethn==''):
        field=med+'_'+'AC'+'_'+sex
    elif (sex==''):
        field=med+'_'+'AC'+'_'+ethn
    else:
        field=med+'_'+'AC'+'_'+ethn+'_'+sex
    
    if valueOfField(field,string) == -1:
        return False
    else:
        return True
    
def popName(med, ethn, sex):
    if (med=='') & (ethn==''):
        field=sex
    elif (ethn=='') & (sex==''):
        field=med
    elif (med=='') & (sex==''):
        field=ethn
    elif (med==''):
        field=ethn+'_'+sex
    elif (ethn==''):
        field=med+'_'+sex
    elif (sex==''):
        field=med+'_'+ethn
    else:
        field=med+'_'+ethn+'_'+sex
        
    return field

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



def populationInfo(med, quant, ethn, sex, string):
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
        
    return valueOfField(field,string)


def listOfAvailablePopulations(string):
    medList = ['controls', 'non_cancer', 'non_neuro', 'non_topmed', '']
    ethnList = ['afr', 'amr', 'nfe', 'sas', 'fin', 'asj', 'eas', 'nfe_seu', 'eas_kor', 'nfe_bgr', 'nfe_swe',
                 'eas_jpn', 'nfe_nwe', 'nfe_onf', 'nfe_est', 'oth', '']
    sexList=['female', 'male', '']
    listPop =[]
    
    for i in medList:
        for j in ethnList:
            for k in sexList:
                if populationPresent(i,j,k,string)==True:
                    listPop.append(popName(i,j,k))
    return listPop
   
        
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def listToCSV(listName, csvName):
    rows=[[i] for i in listName]
    with open(csvName, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)

def csvLinesToListOfFloats(csvName):
    csvFile= open(csvName,'r')
    csvLines= csvFile.readlines()

    listFloats=[]
    for i in range(len(csvLines)):
        listFloats.append(float(csvLines[i][:-1]))

    return listFloats

def VCFtoCSV(vcf_file,infoList,nbInd,csv_file):

    header=[]
    
    header.append('chromo')
    header.append('position')
    header.append('ID')
    header.append('ref')
    header.append('alt')
    header.append('qual')
    header.append('filter')

    for field in infoList:
        header.append(field)

    header.append('AF')

    header.append('Observed_Hete')
    header.append('Expected_Hete')
    header.append('Excess_Hete')
    
    header.append('Observed_HomAlt')
    header.append('Expected_HomAlt')
    header.append('Excess_HomAlt')

    header.append('Observed_HomRef')
    header.append('Expected_HomRef')
    header.append('Excess_HomRef')
    
    header.append('pv')
    
    data=[]
    
    vcfFile= open(vcf_file,'r')
    vcfLines= vcfFile.readlines()

    for l in vcfLines:    

        if l[0] != '#':
            #print(rr)
            dataLine=[]

            lineSplit=l.split()

            for i in range(7):
                dataLine.append(lineSplit[i])
            
            info=lineSplit[7]
            
            for field in infoList:
                dataLine.append(valueOfField(field,info))

            n=nbInd
            
            nA = int(valueOfField('AC',info))

            dataLine.append(nA/(2*n))
            
            #-------------------------Expected vs Observed----------------------------
            genotypes = lineSplit[9:len(lineSplit)]
            
            nAA = genotypes.count('1|1')
            nAB = nA-2*nAA  
            nBB = n - nAA - nAB

            Expected_nAA = ExpecteHomAlt(n,nA)
            Expected_nAB = ExpectedHete(n,nA)
            Expected_nBB = n - Expected_nAA - Expected_nAB
            

            dataLine.append(nAB)
            dataLine.append(Expected_nAB)
            if nAB > Expected_nAB:
                dataLine.append(1)
            else:
                dataLine.append(0)                        
                
            dataLine.append(nAA)
            dataLine.append(Expected_nAA)
            if nAA > Expected_nAA:
                dataLine.append(1)
            else:
                dataLine.append(0)            

            dataLine.append(nBB)
            dataLine.append(Expected_nBB)
            if nBB > Expected_nBB:
                dataLine.append(1)
            else:
                dataLine.append(0)


            
            #-------------------------p-value----------------------------
            
            if nA==0:
                dataLine.append('error')
            else:
                dataLine.append(midPvalue(LH_Distrib(n,nA),nAB))
            
            data.append(dataLine)

            
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)        


def Distrib_list(inputList,x,lowerThreshold,col):   
              
    NbPerBags=[0]*len(x)
    NbPerBagsToPrint=[0]*len(x)
    
    for i in inputList:    
        for j in range(0,len(x)):
        
            if j==0:
                if lowerThreshold<i and i<=x[0]:
                    NbPerBags[j]+=1
        
            else:
                if x[j-1]<i and i<=x[j]:
                    NbPerBags[j]+=1
                            
    print("total number of variants plotted is "+ str(sum(NbPerBags)));print('\n')
    total=sum(NbPerBags)
    print(NbPerBags);print('\n')
    for i in range(len(NbPerBags)):
        NbPerBags[i]=NbPerBags[i]/total*100
        NbPerBagsToPrint[i]=round(NbPerBags[i],2)
    print(NbPerBagsToPrint);print('\n')
    
    
    figure(figsize=(4, 3), dpi=120)
    
    x_str=[]
    for i in range(len(x)):
        if i == 0:
            # x_str.append('< 0 ')
            x_str.append('0 - '+str(x[i]))
        else:
            x_str.append(str(x[i-1]) + ' - '+str(x[i]))
    
    #print(NbPerBags)
    plt.ylabel("% of variants  ",fontsize='10')  
    plt.xticks(fontsize='8',rotation=45)
    plt.yticks(fontsize='8')
    #plt.title("Disitribution of variants' MAF in "+ pop +' population with MAF above ' + str(lowerThreshold))
    plt.bar(x_str, NbPerBags, color=col)
    plt.show()


def statsInfo(csvfile,nbInd,extra):
    retrieve(csvfile,nbInd);
    
    obsAlleFreq=csvLinesToListOfFloats(csvfile[:-4]+'_obsAlleFreq.csv')
    obsHeteFreq=csvLinesToListOfFloats(csvfile[:-4]+'_obsHeteFreq.csv')
    expHeteFreq=csvLinesToListOfFloats(csvfile[:-4]+'_expHeteFreq.csv')
    PV_HWE=csvLinesToListOfFloats(csvfile[:-4]+'_PV_HWE.csv')
    PO=csvLinesToListOfFloats(csvfile[:-4]+'_PO.csv')


    ######----------------------------------------------------------------------------------------
    
    indexAncest =[]
    indexDeriv =[]

    obsAlleFreq_Ancest=[]
    obsAlleFreq_Derived=[]

    obsHeteFreq_Ancest=[]
    obsHeteFreq_Derived=[]

    PV_HWE_Ancest=[]
    PV_HWE_Derived=[]


    for i in range(len(PO)):
        if int(PO[i]) == -1:
            indexAncest.append(i)
            obsAlleFreq_Ancest.append(obsAlleFreq[i])
            obsHeteFreq_Ancest.append(obsHeteFreq[i])
            PV_HWE_Ancest.append(PV_HWE[i])
        else:
            indexDeriv.append(i)
            obsAlleFreq_Derived.append(obsAlleFreq[i])
            obsHeteFreq_Derived.append(obsHeteFreq[i]) 
            PV_HWE_Derived.append(PV_HWE[i]) 

    t=np.arange(0, 1, 0.01)
    t_1=np.arange(0.05, 0.5, 0.01)
    t_2=np.arange(0.5, 0.95, 0.01)

    ######----------------------------------------All variants + HWE ----------------------------------------------------
    print('------------------------------------------------------------------------------------------------');print('\n')
    print('All variants + HWE');print('\n')
    print('Nb of variants: ' + str(len(obsAlleFreq)));print('\n')

    print('AF distribution');print('\n')
    Distrib_list(obsAlleFreq,[0.001,0.01,0.1,1],-0.0001,'blue')

    print('p-value distribution');print('\n')
    Distrib_list(PV_HWE,[0.001,0.01,0.1,1],-0.0001,'blue')
    PV_inflation(PV_HWE,[0.001,0.01,0.05])

    
    plt.figure(figsize=(11, 7), dpi=120)
    plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
    plt.hexbin(obsAlleFreq, obsHeteFreq, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
    plt.colorbar()
    plt.show()

    computingAlpha(obsAlleFreq,obsHeteFreq)

    if extra == 'yes':
        ######----------------------------------------Ancestral variants ----------------------------------------------------
        print('------------------------------------------------------------------------------------------------');print('\n')
        print('Ancestral variants');print('\n')
        print('Nb of variants: ' + str(len(obsAlleFreq_Ancest)));print('\n')
    
        print('AF distribution');print('\n')
        Distrib_list(obsAlleFreq_Ancest,[0.001,0.01,0.1,1],-0.0001,'blue')

        print('p-value distribution');print('\n')
        Distrib_list(PV_HWE_Ancest,[0.001,0.01,0.1,1],-0.0001,'blue')
        PV_inflation(PV_HWE_Ancest,[0.001,0.01,0.05])
        
        plt.figure(figsize=(11, 7), dpi=120)
        plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
        plt.hexbin(obsAlleFreq_Ancest, obsHeteFreq_Ancest, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
        plt.colorbar()
        plt.show()

        computingAlpha(obsAlleFreq_Ancest,obsHeteFreq_Ancest)
    
    
        ######----------------------------------------Derived variants ----------------------------------------------------
        print('------------------------------------------------------------------------------------------------');print('\n')
        print('Derived variants');print('\n')
        print('Nb of variants: ' + str(len(obsAlleFreq_Derived)));print('\n')

        print('AF distribution');print('\n')
        Distrib_list(obsAlleFreq_Derived,[0.001,0.01,0.1,1],-0.0001,'blue')

        print('p-value distribution');print('\n')
        Distrib_list(PV_HWE_Derived,[0.001,0.01,0.1,1],-0.0001,'blue')
        PV_inflation(PV_HWE_Derived,[0.001,0.01,0.05])
        
        plt.figure(figsize=(11, 7), dpi=120)
        plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
        plt.hexbin(obsAlleFreq_Derived, obsHeteFreq_Derived, norm=mpl.colors.LogNorm(), gridsize=(15,15), cmap=plt.cm.jet, mincnt=1)
        plt.colorbar()
        plt.show()

        computingAlpha(obsAlleFreq_Derived,obsHeteFreq_Derived)

def retrieve(variantsCSV,nbInd):
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()

    obsAlleFreq = []
    obsHeteFreq = []
    expHeteFreq = []
    listPV_HWE = []
    list_PO = []

    for i in range(1,len(csvLines)):

        info = csvLines[i][:-1].split(",")


        n=nbInd
        nAB=int(info[16])
        AF=float(info[15])

        obsAlleFreq.append(AF)
        obsHeteFreq.append(nAB/n)
        expHeteFreq.append(2*AF*(1-AF))
        listPV_HWE.append(float(info[25]))
        list_PO.append(info[10])


    listToCSV(obsAlleFreq, variantsCSV[:-4]+'_obsAlleFreq.csv')
    listToCSV(obsHeteFreq, variantsCSV[:-4]+'_obsHeteFreq.csv')
    listToCSV(expHeteFreq, variantsCSV[:-4]+'_expHeteFreq.csv')
    listToCSV(listPV_HWE, variantsCSV[:-4]+'_PV_HWE.csv')
    listToCSV(list_PO, variantsCSV[:-4]+'_PO.csv')


def computingAlpha(obsAlleFreq,obsHeteFreq):
    small_x = np.array(obsAlleFreq) 
    Y = np.array(obsHeteFreq)
    X = np.subtract(np.power(small_x, 2), small_x)
    alpha = np.dot(X,X)**(-1) * np.dot(X,Y)

    print('Nb of variants is ' + str(len(obsAlleFreq)))
    print('\n')
    print("the parameter alpha is equal to " + str(alpha))
    print('\n')

    return alpha

def computingAlpha_2(obsAlleFreq,obsHeteFreq,PV_HWE,tail):

    obsAlleFreq_new=[]
    obsHeteFreq_new=[]

    lower=np.percentile(PV_HWE, tail)
    
    for l in range(len(obsAlleFreq)):
        if PV_HWE[l] > lower:
            obsAlleFreq_new.append(obsAlleFreq[l])
            obsHeteFreq_new.append(obsHeteFreq[l]) 
    
    small_x = np.array(obsAlleFreq_new) 
    Y = np.array(obsHeteFreq_new)
    X = np.subtract(np.power(small_x, 2), small_x)
    alpha = np.dot(X,X)**(-1) * np.dot(X,Y)


    print('Nb of variants is ' + str(len(obsAlleFreq)))
    print('\n')
    print('Nb of variants who passed the QC is ' + str(len(obsAlleFreq_new)))
    print('\n')
    print("the parameter alpha is equal to " + str(alpha))
    print('\n')

    return alpha


def statsInfo_2(csvfile,nbInd,extra,tail):
    retrieve(csvfile,nbInd);
    
    obsAlleFreq=csvLinesToListOfFloats(csvfile[:-4]+'_obsAlleFreq.csv')
    obsHeteFreq=csvLinesToListOfFloats(csvfile[:-4]+'_obsHeteFreq.csv')
    expHeteFreq=csvLinesToListOfFloats(csvfile[:-4]+'_expHeteFreq.csv')
    PV_HWE=csvLinesToListOfFloats(csvfile[:-4]+'_PV_HWE.csv')
    PO=csvLinesToListOfFloats(csvfile[:-4]+'_PO.csv')


    ######----------------------------------------------------------------------------------------
    
    indexAncest =[]
    indexDeriv =[]

    obsAlleFreq_Ancest=[]
    obsAlleFreq_Derived=[]

    obsHeteFreq_Ancest=[]
    obsHeteFreq_Derived=[]

    PV_HWE_Ancest=[]
    PV_HWE_Derived=[]


    for i in range(len(PO)):
        if int(PO[i]) == -1:
            indexAncest.append(i)
            obsAlleFreq_Ancest.append(obsAlleFreq[i])
            obsHeteFreq_Ancest.append(obsHeteFreq[i])
            PV_HWE_Ancest.append(PV_HWE[i])
        else:
            indexDeriv.append(i)
            obsAlleFreq_Derived.append(obsAlleFreq[i])
            obsHeteFreq_Derived.append(obsHeteFreq[i]) 
            PV_HWE_Derived.append(PV_HWE[i]) 


    ######----------------------------------------------------------------------------------------
    
    indexHWD =[]
    indexHWE =[]

    obsAlleFreq_HWD=[]
    obsAlleFreq_HWE=[]

    obsHeteFreq_HWD=[]
    obsHeteFreq_HWE=[]


    for i in range(len(PV_HWE)):
        if PV_HWE[i] < 1e-2:
            indexHWD.append(i)
            obsAlleFreq_HWD.append(obsAlleFreq[i])
            obsHeteFreq_HWD.append(obsHeteFreq[i])
        else:
            indexHWE.append(i)
            obsAlleFreq_HWE.append(obsAlleFreq[i])
            obsHeteFreq_HWE.append(obsHeteFreq[i])

    ######------------------------------------------------------------------------------------------ 
    obsAlleFreq_new=[]
    obsHeteFreq_new=[]

    lower=np.percentile(PV_HWE, 2.5)
    
    for l in range(len(obsAlleFreq)):
        if PV_HWE[l] > lower:
            obsAlleFreq_new.append(obsAlleFreq[l])
            obsHeteFreq_new.append(obsHeteFreq[l])

    ######------------------------------------------------------------------------------------------

    t=np.arange(0, 1, 0.01)
    t_1=np.arange(0.05, 0.5, 0.01)
    t_2=np.arange(0.5, 0.95, 0.01)

    ######----------------------------------------All variants + HWE ----------------------------------------------------
    print('------------------------------------------------------------------------------------------------');print('\n')
    print('All variants + HWE');print('\n')
    print('Nb of variants: ' + str(len(obsAlleFreq)));print('\n')

    print('AF distribution');print('\n')
    Distrib_list(obsAlleFreq,[0,0.001,0.01,0.1,1],-0.0001,'blue')

    print('p-value distribution');print('\n')
    Distrib_list(PV_HWE,[0,0.001,0.01,0.1,1],-0.0001,'blue')
    PV_inflation(PV_HWE,[0.001,0.01,0.05])

    plt.figure(figsize=(11, 7), dpi=120)
    plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
    plt.hexbin(obsAlleFreq, obsHeteFreq, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
    plt.colorbar()
    plt.show()
    

    ######---------------------------------------- Variants in HWD + HWE ---------------------------------------------------
    print('------------------------------------------------------------------------------------------------');print('\n')
    print('Variants in HWD + HWE');print('\n')
    print('Nb of variants: ' + str(len(obsAlleFreq_HWD)));print('\n')
    
    
    plt.figure(figsize=(11, 7), dpi=120)
    plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
    plt.hexbin(obsAlleFreq_HWD, obsHeteFreq_HWD, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
    plt.colorbar()
    plt.show()

    ######---------------------------------------- All variants + HWE + CE ---------------------------------------------------
    print('------------------------------------------------------------------------------------------------');print('\n')
    print('All variants + HWE + CE');print('\n')
    print('Nb of variants: ' + str(len(obsAlleFreq)));print('\n')

    alpha = -computingAlpha(obsAlleFreq,obsHeteFreq)
    
    plt.figure(figsize=(11, 7), dpi=120)
    plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
    plt.plot(t, alpha*t*(1-t), '--', color = 'lime', alpha = 1)
    plt.hexbin(obsAlleFreq, obsHeteFreq, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
    plt.colorbar()
    plt.show()


    addCE(csvfile,alpha,csvfile[:-4]+'_MORE.csv',nbInd)
    retrieve_CE(csvfile[:-4]+'_MORE.csv')
    statsInfo_CE(csvfile,nbInd,'no',alpha,2.5)

    


#     ######----------------------------------------Variants that passed QC + HWE + CE ----------------------------------------------------
#     print('------------------------------------------------------------------------------------------------');print('\n')
#     print('All variants that passed QC + HWE + CE');print('\n')
#     print('Nb of variants: ' + str(len(obsAlleFreq)));print('\n')


#     alpha = -computingAlpha_2(obsAlleFreq,obsHeteFreq,PV_HWE,tail)

#     plt.figure(figsize=(11, 7), dpi=120)
#     plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
#     plt.plot(t, alpha*t*(1-t), '--', color = 'lime', alpha = 1)
#     plt.hexbin(obsAlleFreq_new, obsHeteFreq_new, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
#     plt.colorbar()
#     plt.show()

#     addCE(csvfile,alpha,csvfile[:-4]+'_MORE.csv',nbInd)
#     retrieve_CE(csvfile[:-4]+'_MORE.csv')
#     statsInfo_CE(csvfile,nbInd,'no',alpha,2.5)
    

    if extra == 'yes':
        ######----------------------------------------Ancestral variants ----------------------------------------------------
        print('------------------------------------------------------------------------------------------------');print('\n')
        print('Ancestral variants');print('\n')
        print('Nb of variants: ' + str(len(obsAlleFreq_Ancest)));print('\n')
    
        print('AF distribution');print('\n')
        Distrib_list(obsAlleFreq_Ancest,[0,0.001,0.01,0.1,1],-0.0001,'blue')

        print('p-value distribution');print('\n')
        Distrib_list(PV_HWE_Ancest,[0,0.001,0.01,0.1,1],-0.0001,'blue')
        PV_inflation(PV_HWE_Ancest,[0.001,0.01,0.05])
        
        plt.figure(figsize=(11, 7), dpi=120)
        plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
        plt.hexbin(obsAlleFreq_Ancest, obsHeteFreq_Ancest, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
        plt.colorbar()
        plt.show()

        computingAlpha(obsAlleFreq_Ancest,obsHeteFreq_Ancest)
    
    
        ######----------------------------------------Derived variants ----------------------------------------------------
        print('------------------------------------------------------------------------------------------------');print('\n')
        print('Derived variants');print('\n')
        print('Nb of variants: ' + str(len(obsAlleFreq_Derived)));print('\n')

        print('AF distribution');print('\n')
        Distrib_list(obsAlleFreq_Derived,[0,0.001,0.01,0.1,1],-0.0001,'blue')

        print('p-value distribution');print('\n')
        Distrib_list(PV_HWE_Derived,[0,0.001,0.01,0.1,1],-0.0001,'blue')
        PV_inflation(PV_HWE_Derived,[0.001,0.01,0.05])
        
        plt.figure(figsize=(11, 7), dpi=120)
        plt.plot(t, 2*t*(1-t), '--', color = 'red', alpha = 1)
        plt.hexbin(obsAlleFreq_Derived, obsHeteFreq_Derived, norm=mpl.colors.LogNorm(), gridsize=(15,15), cmap=plt.cm.jet, mincnt=1)
        plt.colorbar()
        plt.show()

        computingAlpha(obsAlleFreq_Derived,obsHeteFreq_Derived)


def addCE(variantsCSV,alpha,newCSV,nbInd):

    data = []
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()

    header=csvLines[0][:-1].split(",")
    header.append('PV_CE')

    header.append('Expected_Hete_CE')
    header.append('Excess_Hete_CE')

    header.append('Expected_HomAlt_CE')
    header.append('Excess_HomAlt_CE')

    header.append('Expected_HomRef_CE')
    header.append('Excess_HomRef_CE')

    header.append('Prob_HWE')
    header.append('Prob_CE')

    data.append(header)


    for i in range(1,len(csvLines)):
        #if i%100000==0:
        #    print(str(i) + ' lines read at ' + str(datetime.now()))

        info = csvLines[i][:-1].split(",")

        nA=int(info[13])
        n=nbInd
        nAB = int(info[16])
        pA = float(nA/(2*n))
        theta = degreeOfDiseq(alpha,pA)

        #----------------------------------PV_CE --------------------

        PV_CE = midPvalue(CE_Distrib(n,nA,theta,alpha),nAB)



        info.append(PV_CE)

        #----------------------------------CE expected  --------------------

        obsHomo=int(info[19])
        obsHete = nbHete(nA,obsHomo)
        obsHomRef=n-obsHomo-obsHete

        expHomo = ExpecteHomAlt_CE(n,nA, alpha)
        expHete = ExpectedHete_CE(n,nA, alpha)
        expHomoRef = n - expHete - expHomo


        info.append(expHete)
        if obsHete > expHete:
            info.append(1)
        else:
            info.append(0)

        info.append(expHomo)
        if obsHomo > expHomo:
            info.append(1)
        else:
            info.append(0)

        info.append(expHomoRef)
        if obsHomRef > expHomoRef:
            info.append(1)
        else:
            info.append(0)


        #----------------------------------Probs under HWE and CE  --------------------
        Prob_HWE = LH_Distrib(n,nA)[nAB]
        Prob_CE = CE_Distrib(n,nA,theta,alpha)[nAB]

        info.append(Prob_HWE)
        info.append(Prob_CE)

        data.append(info)

    #print("\n");print("STARTS WRITING CSV", datetime.now());print("\n")
    with open(newCSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)


def retrieve_CE(variantsCSV):
    csvFile= open(variantsCSV,'r')
    csvLines= csvFile.readlines()

    listPV_CE = []

    for i in range(1,len(csvLines)):

        info = csvLines[i][:-1].split(",")
        listPV_CE.append(float(info[26]))

    listToCSV(listPV_CE, variantsCSV[:-4]+'_PV_CE.csv')


def statsInfo_CE(csvfile,nbInd,extra,alpha,tail):
    retrieve(csvfile,nbInd);
    
    obsAlleFreq=csvLinesToListOfFloats(csvfile[:-4]+'_obsAlleFreq.csv')
    obsHeteFreq=csvLinesToListOfFloats(csvfile[:-4]+'_obsHeteFreq.csv')
    expHeteFreq=csvLinesToListOfFloats(csvfile[:-4]+'_expHeteFreq.csv')
    PV_HWE=csvLinesToListOfFloats(csvfile[:-4]+'_PV_HWE.csv')
    PO=csvLinesToListOfFloats(csvfile[:-4]+'_PO.csv')
    PV_CE=csvLinesToListOfFloats(csvfile[:-4]+'_MORE_PV_CE.csv')


    ######------------------------------------------------------------------------------------------       
            
    indexCD =[]
    indexCE =[]

    obsAlleFreq_CD=[]
    obsAlleFreq_CE=[]

    obsHeteFreq_CD=[]
    obsHeteFreq_CE=[]


    for i in range(len(PV_CE)):
        if PV_CE[i] < 1e-2:
            indexCD.append(i)
            obsAlleFreq_CD.append(obsAlleFreq[i])
            obsHeteFreq_CD.append(obsHeteFreq[i])
        else:
            indexCE.append(i)
            obsAlleFreq_CE.append(obsAlleFreq[i])
            obsHeteFreq_CE.append(obsHeteFreq[i])


    ######------------------------------------------------------------------------------------------ 
    obsAlleFreq_new=[]
    obsHeteFreq_new=[]

    lower=np.percentile(PV_HWE, 2.5)
    
    for l in range(len(obsAlleFreq)):
        if PV_HWE[l] > lower:
            obsAlleFreq_new.append(obsAlleFreq[l])
            obsHeteFreq_new.append(obsHeteFreq[l])

    ######----------------------------------------------------------------------------------------

    t=np.arange(0, 1, 0.01)
    t_1=np.arange(0.05, 0.5, 0.01)
    t_2=np.arange(0.5, 0.95, 0.01)

    ######----------------------------------------All variants + CE ----------------------------------------------------
    print('------------------------------------------------------------------------------------------------');print('\n')
    print('All variants + CE');print('\n')
    print('Nb of variants: ' + str(len(obsAlleFreq)));print('\n')
    print("The parameter alpha is equal to " + str(alpha))
    print('\n')


    print('p-value distribution');print('\n')
    Distrib_list(PV_CE,[0,0.001,0.01,0.1,1],-0.0001,'blue')
    PV_inflation(PV_CE,[0.001,0.01,0.05])

    plt.figure(figsize=(11, 7), dpi=120)
    plt.plot(t, (alpha)*t*(1-t), '--', color = 'lime', alpha = 1)
    plt.hexbin(obsAlleFreq, obsHeteFreq, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
    plt.colorbar()
    plt.show()



    ######---------------------------------------- Variants in CD + CE ---------------------------------------------------
    print('------------------------------------------------------------------------------------------------');print('\n')
    print('Variants in CD + CE');print('\n')
    print('Nb of variants: ' + str(len(obsAlleFreq_CD)));print('\n')
    
    
    plt.figure(figsize=(11, 7), dpi=120)
    plt.plot(t, alpha*t*(1-t), '--', color = 'lime', alpha = 1)
    plt.hexbin(obsAlleFreq_CD, obsHeteFreq_CD, norm=mpl.colors.LogNorm(), gridsize=(75,75), cmap=plt.cm.jet, mincnt=1)
    plt.colorbar()
    plt.show()


def PV_inflation(inputList,x):   
                            
    print("total number of variants  is "+ str(len(inputList)));print('\n')

    for thr in x:
        nbVar = len([i for i in inputList if i <= thr])
        percVar = round(nbVar*100/len(inputList),2)
        print("Nb of variants with p-values less than " + str(thr*100) + "%: " + str(nbVar) + " (" + str(percVar) +"%)");print('\n')


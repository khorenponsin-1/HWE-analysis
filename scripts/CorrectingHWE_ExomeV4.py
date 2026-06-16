import pandas as pd
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import os.path
from os import path
import math
import numpy as np

from HWE import midPvalue, ExpectedHete, CE_Distrib, degreeOfDiseq, nbHete, ExpecteHomAlt_CE, ExpectedHete_CE, LH_Distrib


def addCE_ExomeV4(variantsCSV,alpha,newCSV):

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

    header.append('H_o/H_e')
    header.append('(H_o+1)/(H_e+1)')

    data.append(header)

    for line in csvLines[1:len(csvLines)]:
#     for i in range(1,len(csvLines)):
        #if i%100000==0:
        #    print(str(i) + ' lines read at ' + str(datetime.now()))

        info = line[:-1].split(",")

        nA=int(info[8])
        n=int(info[9])/2
        nAB = int(info[113])
        pA = float(nA/(2*n))
        theta = degreeOfDiseq(alpha,pA)

        #----------------------------------PV_CE --------------------

        PV_CE = midPvalue(CE_Distrib(n,nA,theta,alpha),nAB)



        info.append(PV_CE)

        #----------------------------------CE expected  --------------------

        obsHomo=int(info[11])
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


        #----------------------------------Add Observed / Expected  --------------------
        ratio = (obsHete)/(expHete)
        ratio_plus1 = (obsHete+1)/(expHete+1)

        info.append(ratio)
        info.append(ratio_plus1)




        data.append(info)

    #print("\n");print("STARTS WRITING CSV", datetime.now());print("\n")
    with open(newCSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)

    csvFile.close()
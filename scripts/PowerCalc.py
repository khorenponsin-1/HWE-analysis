import pandas as pd
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import os.path
from os import path
import math
import numpy as np
import plotly.express as px


from HWE import midPvalue, ExpectedHete, CE_Distrib, degreeOfDiseq, nbHete, ExpecteHomAlt_CE, ExpectedHete_CE, LH_Distrib


def degreeOfDiseq_2(pAA,pAB,pBB):
    if pAA == 0 or pBB == 0:
        return 1e100
    elif pAB == 0:
        return 1e-100
    else:
        return pAB**2/(pAA*pBB)
    
    
def typeIError_midPvalue(probs,alpha):
    
    typeIError = 0
    
    for i in range(len(probs)):
        if probs[i] > 1e-10:
            if midPvalue(probs,i) < alpha:
                typeIError += probs[i]
    
    return typeIError

def power_test(probs_H0,probs_H1,alpha):
    
    power = 0
    
    for i in range(len(probs_H0)):
        if probs_H1[i] > 1e-10:
            if midPvalue(probs_H0,i) < alpha:
                power += probs_H1[i]
    
    return power


def CE_Distrib_2(n,nA,theta,nAB):
    
    nB=2*n-nA
    expected = nAB
    
    
    probs=[0]*(nA+1)
    scale=1
    probs[expected]=scale
            
    #Recurrence all the way down
    nAB=expected
    nAA=int((nA-nAB)/2)
    nBB=n-nAB-nAA
    for i in range(expected,1,-2):
        
        probs[i-2]=probs[i]*nAB*(nAB-1)/(theta*(nAA+1)*(nBB+1))
        nAB-=2
        nAA+=1
        nBB+=1
          
    #Recurrence all the way up
    nAB=expected
    nAA=int((nA-nAB)/2)
    nBB=n-nAB-nAA
    for i in range(expected, nA-1,2):
        
        probs[i+2]=probs[i]*theta*nAA*nBB/((nAB+2)*(nAB+1))
        nAB+=2
        nAA-=1
        nBB-=1       
    
    #Scaling to retrieve probability values
    totalCoef=sum(probs)
    for i in range(len(probs)):
        probs[i] = probs[i]/totalCoef
         
    return(probs)


def test(pA,n,sAA,sAB,sBB,toPrint,alpha):
    pAA = pA**2
    pAB = 2*pA*(1-pA)
    pBB = (1-pA)**2
    
#     print(pAB**2/(pAA*pBB))
    
    pAA_abs = sAA * pAA
    pAB_abs = sAB * pAB
    pBB_abs = sBB * pBB
    
    p_abs_tot = pAA_abs + pAB_abs + pBB_abs
    
    pAA_new = pAA_abs / p_abs_tot
    pAB_new = pAB_abs / p_abs_tot
    pBB_new = pBB_abs / p_abs_tot
    
    pA_new = pAA_new + pAB_new/2
        
    nAA_Obs = round(pAA_new * n)
    nAB_Obs = round(pAB_new * n)
    nBB_Obs = round(pBB_new * n)
    
    nA = nAA_Obs*2 + nAB_Obs
    
    pv_HWE = midPvalue(LH_Distrib(n,nA),nAB_Obs)
    
    theta = degreeOfDiseq_2(pAA_new,pAB_new,pBB_new)
    
    H0 = LH_Distrib(n,nA)
    H1 = CE_Distrib_2(n,nA,theta,nAB_Obs)

    
    power = power_test(H0,H1,alpha)
    
    if toPrint == 'yes':
        print('#--------------------- Pre sample selection AF = ' + str(round(pA,4)) + ' ---------------------');print('\n')
#         print('Pre selection AF = ' + str(round(pA,4)));print('\t')
#         print('Post selection AF = ' + str(round(pA_new,4)));print('\t')
        print('sAA = ' + str(sAA));print('\t')
        print('sAB = ' + str(sAB));print('\t')
        print('sBB = ' + str(sBB));print('\t')
        print('Pre selection AF = ' + str(pA));print('\t')
        print('Post selection AF = ' + str(pA_new));print('\t')
        print('With ' + str(n) + ' individuals, the p-value in post sample selection dataset is ' + str(pv_HWE));print('\t')
        print('Power is equal to ' + str(power));print('\t')
        print('Theta is equal to ' + str(theta));print('\t')
        
        if len(H0)%2 == 0:
            H0_toPlot = H0[1::2]
            H0_NbHeteToPlot = np.arange(len(H0))[1::2]
            H1_toPlot = H1[1::2]
            H1_NbHeteToPlot = np.arange(len(H1))[1::2]
        else:
            H0_toPlot = H0[::2]
            H0_NbHeteToPlot = np.arange(len(H0))[::2]
            H1_toPlot = H1[::2]
            H1_NbHeteToPlot = np.arange(len(H1))[::2]
            
        plt.figure(figsize=(20,10))
        plt.plot(H0_NbHeteToPlot, H0_toPlot, marker = '.',linestyle='-',color = 'g')
        plt.plot(H1_NbHeteToPlot, H1_toPlot, marker = '.',linestyle='-',color = 'r')
        plt.axhline(y = 0, color = 'black', linestyle = '-',linewidth = 5) 
        plt.show()
        plt.close()
        
    return pA_new, pv_HWE, power




    
    

def violin_plot(values, null, alt, ax):
    
    samples_null = np.random.choice(values, size=100000, p=null)
    samples_alt = np.random.choice(values, size=100000, p=alt)

    violin_parts = ax.violinplot([samples_null, samples_alt], 
                                 showmeans=False, showmedians=False, showextrema=False)

    colors = ['green', 'red']
    for i, violin in enumerate(violin_parts['bodies']):
        violin.set_edgecolor('black')
        violin.set_linewidth(1.5)
        violin.set_facecolor(colors[i])

    ax.set_xticks([1, 2])
    ax.set_xticklabels(['HWE', 'Selection'], fontsize=12)
    
    
def compareGeno(pA, sAA, sAB, sBB, n, asp):

    pAA = pA**2
    pAB = 2*pA*(1-pA)
    pBB = (1-pA)**2

    pAA_abs = sAA * pAA
    pAB_abs = sAB * pAB
    pBB_abs = sBB * pBB

    p_abs_tot = pAA_abs + pAB_abs + pBB_abs

    pAA_new = pAA_abs / p_abs_tot
    pAB_new = pAB_abs / p_abs_tot
    pBB_new = pBB_abs / p_abs_tot

    pA_new = pAA_new + pAB_new/2

    nAA_Obs = round(pAA_new * n)
    nAB_Obs = round(pAB_new * n)
    nBB_Obs = round(pBB_new * n)

    nA = nAA_Obs*2 + nAB_Obs

    pv_HWE = midPvalue(LH_Distrib(n, nA), nAB_Obs)

    theta = degreeOfDiseq_2(pAA_new, pAB_new, pBB_new)

    neutral = LH_Distrib(n, nA)
    selection = CE_Distrib_2(n, nA, theta, nAB_Obs)

    nB = 2*n - nA

    Hete_values = np.arange(len(neutral))
    HomoA_values = []
    HomoB_values = []

    for nbHete in Hete_values:
        HomoA_values.append((nA - nbHete) / 2)
        HomoB_values.append((nB - nbHete) / 2)

    HomoA_values = np.array(HomoA_values)
    HomoB_values = np.array(HomoB_values)

    Hete_values = Hete_values / n
    HomoA_values = HomoA_values / n
    HomoB_values = HomoB_values / n

    print(pA_new)

    fig, axes = plt.subplots(1, 3, figsize=(9, 6), dpi=350)

#     violin_plot(Hete_values, neutral, selection, axes[0])
#     axes[0].set_title('Heterozygotes', fontsize=16, pad=13)
#     axes[0].set_ylabel('Frequency', fontsize=16, labelpad=15)
#     axes[0].tick_params(axis='both', labelsize=14)

    violin_plot(HomoA_values, neutral, selection, axes[0])
    axes[0].set_title('Homozygotes AA', fontsize=16, pad=13)
    axes[0].tick_params(axis='both', labelsize=14)
    
    violin_plot(Hete_values, neutral, selection, axes[1])
    axes[1].set_title('Heterozygotes', fontsize=16, pad=13)
#     axes[0].set_ylabel('Frequency', fontsize=16, labelpad=15)
    axes[1].tick_params(axis='both', labelsize=14)

    violin_plot(HomoB_values, neutral, selection, axes[2])
    axes[2].set_title('Homozygotes BB', fontsize=16, pad=13)
    axes[2].tick_params(axis='both', labelsize=14)

    plt.tight_layout()
    plt.show()
    
    

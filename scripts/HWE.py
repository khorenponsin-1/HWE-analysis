def ExpectedHete(n,nA):
    nB=2*n-nA
    expected = (nA*nB/(2*n))

    if nA%2==0:
        if round(expected)%2==0:
            expected=int(round(expected))

        elif round(expected) > expected:
            expected=round(expected) - 1

        else:
            expected=round(expected) + 1


    else:
        if round(expected)%2==1:
            expected=int(round(expected))

        elif round(expected) > expected:
            expected=round(expected) - 1

        else:
            expected=round(expected) + 1

    return expected

def ExpecteHomAlt(n,nA):
    expected=(nA - ExpectedHete(n,nA))/2
    return round(expected)

def nbHete(nA,altHomo):
    nAB=nA-2*altHomo
    return nAB



def LH_Distrib(n,nA):

    nB=2*n-nA
    expected = ExpectedHete(n,nA)

    probs=[0]*(nA+1)
    scale=1
    probs[expected]=scale

    #Recurrence all the way down
    nAB=expected
    nAA=int((nA-nAB)/2)
    nBB=n-nAB-nAA
    for i in range(expected,1,-2):

        probs[i-2]=probs[i]*nAB*(nAB-1)/(4*(nAA+1)*(nBB+1))
        nAB-=2
        nAA+=1
        nBB+=1

    #Recurrence all the way up
    nAB=expected
    nAA=int((nA-nAB)/2)
    nBB=n-nAB-nAA
    for i in range(expected, nA-1,2):

        probs[i+2]=probs[i]*4*nAA*nBB/((nAB+2)*(nAB+1))
        nAB+=2
        nAA-=1
        nBB-=1

    #Scaling to retrieve probability values
    totalCoef=sum(probs)
    for i in range(len(probs)):
        probs[i] = probs[i]/totalCoef

    return(probs)

def midPvalue(probs,obsAB):


    p_value=0

    probObsAB=probs[obsAB]
    #print(probObsAB)

    for i in range(len(probs)):

        if i == obsAB:
            p_value=p_value + probs[i]/2
            #print('DEDEDEDE')

        elif probs[i] < probObsAB:
            #print(i)
            p_value+=probs[i]

    return p_value

def CE_Distrib(n,nA,theta,alpha):

    nB=2*n-nA
    expected = ExpectedHete_CE(n,nA,alpha)

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

def degreeOfDiseq(alpha,pA):
    pB = 1 - pA
    maf = min(pA,pB)

    pAB = min(2*maf,alpha*pA*pB)

    pAA = pA - pAB/2
    pBB = 1 - pAB - pAA

    if pAA == 0 or pBB == 0:
        return 100
    elif pAB == 0:
        return 0.01
    else:
        return pAB**2/(pAA*pBB)

def ExpectedHete_CE(n,nA, alpha):
    pA = nA/(2*n)
    pB = 1 - pA
    maf = min(pA,pB)

    pAB = min(2*maf,alpha*pA*pB)

    expected = pAB*n

    if nA%2==0:
        if round(expected)%2==0:
            expected=int(round(expected))

        elif round(expected) > expected:
            expected=round(expected) - 1

        else:
            expected=round(expected) + 1


    else:
        if round(expected)%2==1:
            expected=int(round(expected))

        elif round(expected) > expected:
            expected=round(expected) - 1

        else:
            expected=round(expected) + 1

    return expected


def ExpecteHomAlt_CE(n,nA,alpha):
    expected=(nA - ExpectedHete_CE(n,nA,alpha))/2

    return (expected)
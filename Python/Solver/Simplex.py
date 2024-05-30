from numpy import *

def StandardForm(c, boolMin, A, sg, b):
    # Check orders
    n = len(c)
    m = len(b)
    if len(sg) != m or len(A) != m:
        print("ERRO: Matrizes inconsistentes")
        return False, [], []
    for i in range(m):
        if len(A[i]) != n:
            print("ERRO: Matrizes inconsistentes")
            return False, [], []


    # Minimizing
    if boolMin == False:
        c = [(-1)*ci for ci in c]

    # b >= 0
    for i in range(m):
        if b[i] < 0:
            A[i] = [(-1)*a for a in A[i]]
            b[i] = -b[i]
            if sg[i] == ">=":
                sg[i] = "<="
            if sg[i] == "<=":
                sg[i] = ">="

    # Large M
    largeM = 0
    for i in range(m):
        for j in range(n):
            if abs(A[i][j]) > largeM:
                largeM = abs(A[i][j])
    largeM = 1000*largeM
    
    # Slack Variables
    newA = []
    newB = []
    newSg = []
    for i in range(m):
        if sg[i] != "==":
            newA.append(A[i].copy())
            newB.append(b[i])
            newSg.append(sg[i])
        if sg[i] == "==":
            newA.append(A[i].copy())
            newA.append(A[i].copy())
            newB.append(b[i])
            newB.append(b[i])
            newSg.append("<=")
            newSg.append(">=")
    A = newA
    b = newB
    sg = newSg
    m = len(sg)
        
    for i in range(m):
        # Slack variables for <= inequalities
        if sg[i] == "<=":
            n = n+1
            c.append(0)
            for j in range(m):
                if i != j:
                    A[j].append(0)
                else:
                    A[j].append(1)
                sg[i] = "=="
        # Slack variables for >= inequalities
        if sg[i] == ">=":
            n = n+2
            c.append(0)
            c.append(largeM)
            for j in range(m):
                if i != j:
                    A[j].append(0)
                    A[j].append(0)
                else:
                    A[j].append(-1)
                    A[j].append(1)
                sg[i] = "=="



    return c, A, b

# Find Identity Matrix within A
def FindId(A):
    m = len(A)
    n = len(A[0])
    if m > n:
        return False
    
    idMatrix = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        idMatrix[i][i] = 1
    
    indexes = [-1]*n
    for k in range(n):
        column = [A[i][k] for i in range(m)]
        for i in range(m):
            dif = 0
            for j in range(m):
                dif += abs(idMatrix[i][j] - column[j])
            if dif == 0:
                indexes[k] = i
    
    order = [-1]*m
    for i in range(n):
        if indexes[i] > -1:
            order[indexes[i]] = i
    for i in range(m):
        if order[i] == -1:
            return False
    return order

def PrintCurrentTable(c, A, b, bv, cv):
    n = len(c)
    m = len(A)

    print("-"*(22*(n+3)+3))

    print("|", "Basic Variables".center(20, ' '), "|", "Variables".center(20, ' '), end='')
    for i in range(n):
        print("|", f"x[{i}]".center(20, ' '), end='')
    print("|", "Current Value".center(20, ' '), "|")
    
    print("-"*(22*(n+3)+3))

    for j in range(m):
        print("|", f"x[{bv[j]}]".center(20, ' '), "|", f"Constraint {j+1}".center(20, ' '), end='')
        for i in range(n):
            print("|", f"{A[j][i]:0.2f}".center(20, ' '), end='')
        print("|", f"{b[j]:0.2f}".center(20, ' '), "|")

    print("-"*(22*(n+3)+3))

    print("|", " ".center(20, ' '), "|", "Objective Function".center(20, ' '), end='')
    for i in range(n):
        print("|", f"{c[i]:0.2f}".center(20, ' '), end='')
    print("|", f"{cv:0.2f}".center(20, ' '), "|")

    print("-"*(22*(n+3)+3))
    print("")

def InnerProd(u, v):
    n1 = len(u)
    n2 = len(v)
    if (n1 < n2):
        for i in range(n2 - n1):
            u.append(0)
    
    if (n1 > n2):
        for i in range(n1 - n2):
            v.append(0)

    return sum(u[i]*v[i] for i in range(len(u)))


def Simplex(c, boolMin, A, sg, b):
    cSF, ASF, bSF = StandardForm(c, boolMin, A, sg, b)
    n = len(cSF)
    m = len(ASF)
    x = [0]*n
    bv = FindId(ASF)
    obValue = InnerProd(cSF, x)

    PrintCurrentTable(cSF, ASF, bSF, bv, obValue)

    # First check
    canImprove = False
    newBV = -1
    for i in range(n):
        if cSF[i] < 0:
            newBV = i
            canImprove = True
            break
    # Main Simplex Loop
    while canImprove == True:
        # Find BV to replace
        minDiv = 0
        for i in range(m):
            if A[i][newBV] != 0:
                minDiv += abs(b[i]/A[i][newBV])
        minDivLine = -1
        for i in range(m):
            if A[i][newBV] != 0:
                ratio = b[i]/A[i][newBV]
                if ratio > 0 and ratio < minDiv:
                    minDiv = ratio
                    minDivLine = i
        if minDivLine == -1:
            return False, []
        
        # Fixing Constraints in minDivLine
        r = 1/(ASF[minDivLine][newBV])

        bv[minDivLine] = newBV # Setting new basic variable
        bSF[minDivLine] *= r # updating in Current Value
        for j in range(n):
            ASF[minDivLine][j] *= r

        # Fixing Remaining Constraints
        for i in range(m):
            if i != minDivLine:
                currentA = ASF[i][newBV]
                for j in range(n):
                    ASF[i][j] -= ASF[minDivLine][j]*currentA
                bSF[i] -= bSF[minDivLine]*currentA
        
        # Fixing Objective Function
            currentC = cSF[newBV]
            for j in range(n):
                cSF[j] -= ASF[minDivLine][j]*currentC
            obValue -= bSF[minDivLine]*currentC

        PrintCurrentTable(cSF, ASF, bSF, bv, obValue)
        
        # New Verification
        canImprove = False
        newBV = -1
        for i in range(n):
            if cSF[i] < 0:
                newBV = i
                canImprove = True
                break
        

    for i in range(m):
        x[bv[i]] = bSF[i]
    return True, x
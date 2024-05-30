from Solver.Simplex import *

def main():
    filename = "Solver/SimplexInput.txt"
    with open(filename, 'r') as file:
        lines = file.readlines()
        curLine = []

        m = len(lines)-1

        # Função Objetivo
        curLine = lines[0].split()
        n = len(curLine)-1

        c = [0]*n
        A = [[0]*n for j in range(m)]
        sg = [0]*m
        b = [0]*m

        for j in range(n):
            c[j] = float(curLine[j])

        if curLine[n] == 'min':
            boolMin = True
        else:
            boolMin = False

        # Constraints
        for i in range(m):
            curLine = []
            curLine = lines[i+1].split()
            for j in range(n):
                A[i][j] = float(curLine[j])
            sg[i] = curLine[n]
            b[i] = float(curLine[n+1])

    check, x = Simplex(c.copy(), boolMin, A.copy(), sg.copy(), b.copy())

    if check != True:
        print("Não é possível encontrar solução ótima!")
    else:
        print("Solução ótima: {", end='')
        for i in range(len(x)-1):
            print(f"{x[i]:0.2f}", end=', ')
        print(f"{x[len(x)-1]:0.2f}", "}")
        print(f"Valor da Função Objetivo: {InnerProd(c, x):0.2f}")

main()
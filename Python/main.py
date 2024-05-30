from Solver.Simplex import *

def main():
    c = [20, 10]
    boolMin = False
    A = [
        [1,2], [3,1]
    ]
    sg = ["<=", "<="]
    b = [40, 70]

    check, x = Simplex(c.copy(), boolMin, A.copy(), sg.copy(), b.copy())

    if check != True:
        print("Não é possível encontrar solução ótima!")
    else:
        print("Solução ótima:", x)
        print("Valor da Função Objetivo:", InnerProd(c, x))

main()
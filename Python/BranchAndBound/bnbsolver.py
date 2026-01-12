from ortools.linear_solver import pywraplp

class LPSolver:
	def __init__(self, filename):
		with open(filename, 'r') as file:
			lines = file.readlines()
		curLine = []

		# Number of Variables and Constraints
		curLine = lines[0].split()
		self.__n = len(curLine) # Variables
		self.__m = len(lines)-2 # Constraints

		# Initializing matrices
		self.__c = [0]*self.__n
		self.__A = [([0]*self.__n).copy() for j in range(self.__m)]
		self.__b = [0]*self.__m

		# Objective Function Constants
		for j in range(self.__n):
			self.__c[j] = float(curLine[j])

		# Constraints RHS
		curLine = lines[1].split()
		for j in range(self.__m):
			self.__b[j] = float(curLine[j])

		# Constraints
		for i in range(self.__m):
			curLine = []
			curLine = lines[i+2].split()
			for j in range(self.__n):
				self.__A[i][j] = float(curLine[j])

	def GetProblem(self):
		n = 53
		print("-"*n)
		print("[Variáveis]")
		for i in range(self.__n-1):
			print(f"x[{i}], ", end='')
		print(f"x[{self.__n-1}]")
		
		print("[Função Objetivo]")
		print("(min)  ", end='')
		for i in range(self.__n-1):
			print(f"{self.__c[i]}*x[{i}] + ", end='')
		print(f"{self.__c[i]}*x[{self.__n-1}]")
		
		print("[Condições]")
		for i in range(self.__m):
			for j in range(self.__n-1):
				print(f"{self.__A[i][j]}*x[{j}] + ", end='')
			print(f"{self.__A[i][self.__n-1]}*x[{self.__n-1}] <= {self.__b[i]}")
		print("-"*n)

	def RealSolver(self):
		hasSolution, fmax, xRet = self.__GetRealSolver(self.__c, self.__A, self.__b)
		if hasSolution == True:
			print(f"Valor da função objetivo: {fmax:0.3f}")
			print("Solução ótima = [", end="")
			for i in range(len(xRet)-1):
				print(f"{xRet[i]:0.3f}, ", end="")
			print(f"{xRet[len(xRet)-1]:0.3f}]")
		else:
			print("Problema não possui solução ótima.")

	# Real Solver using ORTools
	def __GetRealSolver(self, c, A, b):
		m = len(A)
		n = len(A[0])
		solver = pywraplp.Solver.CreateSolver("GLOP")
		if not solver:
			return

		# [Variables]
		x = [solver.NumVar(0, solver.infinity(), f"x[{i}]") for i in range(n)]


		# [Objective Function]
		solver.Maximize(solver.Sum([c[i]*x[i] for i in range(n)]))
		

		# [Conditions]
		for j in range(m):
			solver.Add(solver.Sum([A[j][i]*x[i] for i in range(n)]) <= b[j])


		status = solver.Solve()
		if status == pywraplp.Solver.OPTIMAL:
			xRet = [x[i].solution_value() for i in range(n)]
			return True, solver.Objective().Value(), xRet
		else:
			return False, [], []
		
	def BnBSolver(self):
		currentNode, bestSolutionValue, bestSolution = self.__BnBSolverRec(self.__c, self.__A, self.__b, 0, -1, [0]*len(self.__c))
		print(f"Número de nós visitados: {currentNode}")
		print(f"Valor da função objetivo: {bestSolutionValue:0.3f}")
		print(f"Solução ótima: {bestSolution}")

	def __BnBSolverRec(self, c, A, b, currentNode, bestSolutionValue, bestSolution):
		sol, fmax, xRet = self.__GetRealSolver(c, A, b)
		m = len(A)
		n = len(A[0])

		print(f"Nó {currentNode}:")

		def NearestInt(v):
			f = int(v)
			if v - f < 0.5:
				return f
			else:
				return f+1

		def CheckIntegerSolution(x):
			approxError = 0.0000001
			for i in range(len(x)):
				if abs(x[i] - NearestInt(x[i])) > approxError:
					return i
			return -1
			
		if sol == False:
			print("Poda por Inviabilidade.")
			return currentNode, bestSolutionValue, bestSolution
		else:
			print(f"Função Objetivo = {fmax:0.3f}")
			print("x = [", end="")
			for i in range(len(xRet)-1):
				print(f"{xRet[i]:0.3f}, ", end="")
			print(f"{xRet[len(xRet)-1]:0.3f}]")

			newCondVar = CheckIntegerSolution(xRet)
			if newCondVar == -1:
				if fmax > bestSolutionValue:
					bestSolutionValue = fmax
					bestSolution = xRet
					print("Poda por Integrabilidade")
				else:
					print("Poda por Optimalidade")
				return currentNode, bestSolutionValue, bestSolution
			else:
				Am = A.copy()
				bm = b.copy()
				line = [0]*n
				line[newCondVar] = 1
				Am.append(line)
				bm.append(int(xRet[newCondVar]))
				currentNode, bestSolutionValue, bestSolution = self.__BnBSolverRec(c, Am, bm, currentNode+1, bestSolutionValue, bestSolution)
				
				AM = A.copy()
				bM = b.copy()
				line[newCondVar] = -1
				AM.append(line)
				bM.append((-1)*(int(xRet[newCondVar])+1))
				currentNode, bestSolutionValue, bestSolution = self.__BnBSolverRec(c, AM, bM, currentNode+1, bestSolutionValue, bestSolution)
				return currentNode, bestSolutionValue, bestSolution
			

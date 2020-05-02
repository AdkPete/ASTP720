import numpy as np
import sys , random
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import scipy.stats as stats
import copy
 
def bin_inc(string , ind):

	'''
	This code takes in a dna string, called string
	Also takes in aan index ind, which should be an index
	for a one or zero in the dna
	This function flips that bit of the string to the other value
	'''

	if ind > len(string):
		print ("error")
		sys.exit()
		
	if ind == 0:
	
		return "1" + string[1::]
	
	elif ind == len(string) - 1:
		return string[0:-1] + "1"
	
	else:
		return string[0:ind] + "1" + string[ind + 1::]

def float_to_bin(N , p):

	'''
	Takes in a float N and precision p
	p should be array - like with length 2
	returns a binary dna string
	has p[0] characters, then a period , then p[1] characters
	'''
	
	if N < 0:
		sign = True
		N *= -1
	else:
		sign = False

	
	i = int(str(N).split(".")[0])
	
	dna_i = '0' * p[0]
	
	for k in range(len(dna_i)):

		if (2 ** (k)) > i:
			break
		elif i % (2 ** (k + 1)) != 0:
			dna_i = bin_inc(dna_i , len(dna_i) - 1 - k)
			i -= 2 ** k
			
	if len(str(N).split(".")) == 0:
		dna_f = '0' * p[1]
		return dna_i + "." + dna_f
	
	f = int(str(N).split(".")[1])
	
	dna_f = "0" * p[1]
	
	
	for k in range(len(dna_f)):

		if (2 ** (k)) > f:
			break
		elif f % (2 ** (k + 1)) != 0:
			dna_f = bin_inc(dna_f , len(dna_f) - 1 - k)
			f -= 2 ** k
	
	
	if not sign:
		
		return "0" + str(dna_i) + "." + dna_f
	
	return "1" + str(dna_i) + "." + dna_f
	
def bin_to_float(N):

	'''
	Takes in a dna string N
	Should be one chromosome
	Converts N to a float
	Returns the float
	'''
	
	if float(N) == 0:
		return 0
	if N[0] == "1":
		sign = -1
	else:
		sign = 1
	
	N = N[1::]
	i = str(N).split(".")[0]
	
	
	a = len(str(i)) - 1
	
	v = 0
	k = 0
	i = str(i)
	
	while a >= 0:
		v += int(i[k]) * (2 ** a)
		a -= 1
		k += 1
		
	if len(str(N).split(".")) == 1:
		return v
	
	f  = str(N).split(".")[1]
	a = len(str(f)) - 1
	
	vf = 0
	k = 0
	f = str(f)
	
	while a >= 0:
		vf += int(f[k]) * (2 ** a)
		a -= 1
		k += 1
		
	res = float(str(v) + "." + str(vf))
	

	return sign * res
	
	
class Genetic:

	def __init__(self , f  , bounds , popsize = 1000, precision = [15 , 15] , creatures = None):
		'''
		f is a function that takes in an array of N parameters and returns some form of goodness of fit / likelihood or similar.
		precision determines how large the dna strings will be
		should be array-like containint 2 elements.
		The first describes the number of values allowed before the decimal point.
		The second describes the number of values allowed after the decimal point.
		'''
		
		self.f = f
		self.prec = precision
		self.bounds = bounds
		self.popsize = popsize
		self.generation = 0
		
		self.best_dna = None
		self.best_f = None
		
		if creatures == None:
			self.creatures = []
		else:
			self.creatures = creatures
			
		
	def create_dna(self , X):
		'''
		X is a numpy array
		Each element should be one of your parameters
		Returns a dna string
		'''
		dna = ""
		for i in X:
			dna += float_to_bin(i , self.prec) + "|"
		
		return dna[0:-1]
		
			
	def initialize(self):
		
		'''
		Sets up a population of creatures for you
		They will be randomly distributed with the specified boundaries
		'''
		
		while len(self.creatures) < self.popsize:
			
			X = np.array( [] )
			for i in self.bounds:
				
				parameter = random.uniform(i[0] , i[1])
				
				X = np.append(X , parameter)
			
			dna = self.create_dna(X)
			self.creatures.append(creature(dna))
			
		
	def determine_fitness(self):
	
		'''
		Measures the fitness for every creature that hasn't had a fitness evaluated yet
		Sorts the list of creatures by fitness
		the first creature will always be the most fit
		'''
		
		for i in self.creatures:
			if i.fitness == None:
			
				c_fit = self.f(i.get_params())
				i.fitness = c_fit
				if self.best_f == None or c_fit > self.best_f:
					self.best_f = c_fit
					self.best_dna = i
				
		self.creatures = sorted(self.creatures , reverse = True)
		
		
	def update(self , N = 1 , makeplots = False):
	
		##Advances our population N times.
		## Defaults to a single update
		## if makeplots = True, this will produce a plot called genumb.pdf for every generation
		## Useful for making movies
		
		for i in range(N):
			
			if makeplots:
				self.nplot()
			if len(self.creatures) < self.popsize:
				self.initialize()
				
			self.determine_fitness()
			self.generation += 1
			
			###PARAMETER
			
			nchildren = int(0.25 * len(self.creatures))
			
			self.creatures = sorted(self.creatures , reverse = True)
			
			child_set = np.array(self.creatures[0:2 * nchildren])
			New_Creatures = []
			while len(child_set) > 0:
				A = random.randint(0 , len(child_set) - 1)
				B = random.randint(0 , len(child_set) - 1)
				while B == A:
					B = random.randint(0 , len(child_set) - 1)
				New = child_set[A].replicate(child_set[B] , self.bounds)
				New.fitness = self.f(New.get_params())
				New.gen = self.generation
				
				New_Creatures.append(New)
				del self.creatures[-1]
				child_set = np.delete(child_set , [A , B] )
				
			for i in New_Creatures:
				self.creatures.append(i)
				
			self.creatures = sorted(self.creatures , reverse = True)
			
			self.best_f = self.creatures[0].fitness
			self.best_dna = self.creatures[0]
		

	def nplot(self):
		
		##Makes a scatter plot of the current distribution of the population
		# Name of save file is the current generation .pdf
		#Produces a pdf, no returns
		
		x0 = []
		x1 = []
		y0 = []
		y1 = []
		for i in self.creatures:
			if i.gen == self.generation:
				x1.append(i.get_params()[0])
				y1.append(i.get_params()[1])
			else:
				x0.append(i.get_params()[0])
				y0.append(i.get_params()[1])
		if self.generation == 0:
			
			plt.scatter(x1 , y1 , color = 'b')
			plt.scatter(x0 , y0 , color = 'b')#color = 'r' , marker = 'x')
		
		else:
			plt.scatter(x0 , y0 , color = 'b')
			plt.scatter(x1 , y1 , color = 'r' , marker = 'x')
		plt.xlim( self.bounds[0][0] , self.bounds[0][1])
		plt.ylim( self.bounds[1][0] , self.bounds[1][1])
		
		plt.savefig(str(self.generation) + ".pdf")
		plt.close()
			
	
		
		
class creature:

	'''
	Class to keep track of useful bits of information for a given creature
	Contains dna, and the functiosn for reproduction / mutation
	Also keeps track of fitness levels
	
	'''
	def __init__(self , dna , gen = 0):
		self.dna = dna
		self.fitness = None
		self.gen = gen
		self.bounds = None
		
		
	def __gt__(self , other):
		if self.fitness > other.fitness:
			return True
		else:
			return False
			
			
	def get_params(self):
	
		'''
		This method just converts this creatures dna into a numpy array filled with the parameters
		'''
		X = np.array( [] )
		for i in self.dna.split("|"):
			X = np.append(X , bin_to_float(i))
		return X
		
	def mutate(self , bounds):
		
		'''
		Mutates the current dna string
		bounds should be a list of lists, telling it which part of the 
		parameter space is reasonable. Will not mutate outside of said bounds
		'''
		
		odna = self.dna
		
		p = 1 / float(len(self.dna))
		
		ndna = ""
		
		for i in range(len(self.dna)):
			if self.dna[i] == "|" or self.dna[i] == ".":
				ndna += self.dna[i]
			else:
				C = np.random.rand()
				if C < p:
					if self.dna[i] == "0":
						ndna += "1"
					else:
						ndna += "0"
				else:
					ndna += self.dna[i]
		self.dna = ndna
		for i in range(len(self.get_params())):
			val = self.get_params()[i]
			if val < bounds[i][0] or val > bounds[i][1]:
				self.dna = odna
				break
			

	def replicate(self , other , bounds):
		
		'''
		Takes in another creature, and a set of bounds for mutations
		Produces and returns a new creature, based on the dna of both parents
		'''
		
		b = random.randint(0 , len(self.dna))

		
		child_dna = ""
		for i in range(len(self.dna)):
		
			if i <= b:
				
			
				child_dna += self.dna[i]
			else:
				
			
				child_dna += other.dna[i]
				

		Result = creature(child_dna)
		Result.mutate(bounds)
		
		return Result

def mkplot(creatures , gens = None):
	
	'''
	Takes in a list of creatures
	Makes a scatter plot of the values of the first two parameters
	'''
	
	if gens == None:
	
		x = []
		y = []
		for i in creatures:
			X = i.get_params()
			x.append(X[0])
			y.append(X[1])
		plt.scatter(x , y)
		plt.xlim(0 , 25)
		plt.ylim(0 , 25)
		plt.show()
		
	else:
		x0 = []
		y0 = []
		
		x1 = []
		y1 = []
		
		for i in creatures:
			if i.gen == gens:
				X = i.get_params()
				x0.append(X[0])
				y0.append(X[1])
			
			else:
				X = i.get_params()
				x1.append(X[0])
				y1.append(X[1])
				
		plt.scatter(x0 , y0 , color = "r" , marker = "x")
		plt.scatter(x1 , y1 , color = "b")
		#plt.xlim(0 , 25)
		#plt.ylim(0 , 25)
		plt.show()
	
def test_fit(X):

	## Sample fitness function

	sig = [ [1.5 , 0 ], [0 , 3.4 ]]
	mu = [ 0 , 0 ]
	
	rv = multivariate_normal( mu , sig )
	
	return rv.pdf(X)	

def mean_fit(X):
	##Harder fitness function, contains a local maximum to avoid
	
	sig = [ [1.5 , 0 ], [0 , 3.4 ]]
	mu = [ 3 , 10 ]
	
	rv = multivariate_normal( mu , sig )
	
	res = rv.pdf(X)
	
	sig = [ [1.5 , 0 ], [0 , 3.4 ]]
	mu = [ 15 , 15 ]
	
	rv = multivariate_normal( mu , sig )
	
	res += 0.5 * rv.pdf(X)
	return res
	
def test_fit_3d(X):
	
	## 3d fitness function
	sig = [ [1.5 , 0 , 0], [0 , 3.4 , 0 ] , [0 , 0 , 2.3] ]
	mu = [ 0 , 0  , 0]

	rv = multivariate_normal( mu , sig )

	return rv.pdf(X)	
	
def mk_subplot(creatures , gen):
	
	'''
	Makes a scatterplot, designed to be used to make multi - paneled plots
	'''
	
	X =  []
	Y =  []
	
	new_X = []
	new_Y = []
	for i in creatures:
		if gen == i.gen:
			new_X.append(i.get_params()[0])
			new_Y.append(i.get_params()[1])
		else:
			X.append(i.get_params()[0])
			Y.append(i.get_params()[1])
			
	plt.scatter(new_X , new_Y , marker = "x" , color = "red")
	plt.scatter(X , Y)
	plt.scatter(0 , 0 , marker = "*" , s = 75 , color = "green")
	
	
def movie():
	'''
	Generates a series of plots that coule be used to make a movie
	'''
	A = Genetic(test_fit_3d , [ [-25 , 25] , [-25 , 25 ] , [-25 , 25] ] , 50 )
	
	for i in range(50):
		
		A.update(makeplots = True)
		
	
	
	

if __name__ == "__main__":

	movie()
	
	#A = Genetic(test_fit_3d , [ [-25 , 25] , [-25 , 25 ] , [-25 , 25] ] , 1000 )
	A = Genetic(test_fit_3d , [ [-25 , 25] , [-25 , 25 ] , [-25 , 25] ] , 150 )
	A.initialize()
	A.determine_fitness()


	#mkplot(A.creatures ,1)
	
	
	for i in range(100):
		if i == 1:
			Gen1 = copy.deepcopy(A.creatures)
		elif i == 5:
			Gen5 = copy.deepcopy(A.creatures)
		
		elif i == 10:
			Gen10 = copy.deepcopy(A.creatures)
		elif i == 15:
			Gen15 = copy.deepcopy(A.creatures)
		
		A.update()

	plt.subplots_adjust(wspace = 0.6 , hspace = 0.6)
	plt.subplot(2 , 2 , 1)
	mk_subplot(Gen1 , 1 )
	plt.title("Generation 1")
	plt.ylabel("y")
	
	plt.subplot(2 , 2 , 2)
	mk_subplot(Gen5 , 5 )
	plt.title("Generation 5")

	
	plt.subplot(2 , 2 , 3)
	mk_subplot(Gen10 , 10 )
	plt.title("Generation 10")
	plt.xlabel("x")
	plt.ylabel("y")
	plt.subplot(2 , 2 , 4)
	mk_subplot(Gen15 , 15 )
	plt.title("Generation 15")
	plt.xlabel("x")
	
	plt.savefig("Test.pdf")
	plt.close()
	
	print (A.creatures[0].get_params())
	#mkplot(A.creatures , 1)

	#print (A.creatures[0].get_params())



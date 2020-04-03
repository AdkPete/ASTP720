import numpy as np
import unittest
import astropy.constants as const
import astropy.units as u
import unittest
import copy , os , shutil , time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


t_ind = 1

def test_mode(nt = 0):
	##Calling this function puts the code into a testing mode
	##Will break some functions. basically assumes that all particles have only one position, not two
	###Useful solely for testing, otherwise do not touch!
	
	global t_ind
	t_ind = nt

def S(r):
	return 1 / np.sqrt(r * r + param.epsilon ** 2)
	
def set_params(fname):
	global param
	param = Params(fname)
	

class Params:

	def __init__(self , paramfile):
		self.paramfile = paramfile
		self.read()
		self.check()
		
	def read(self):
		f = open(self.paramfile)
		for i in f.readlines():
			L = i.split()
			if L[0][0] == "#":
				continue
			elif L[0].lower() == "step_size":
				self.h = float(L[1]) * u.yr
				
			elif L[0].lower() == "end_time":
				self.t_end = float(L[1]) * u.yr
				
			elif L[0].lower() == "epsilon":
				self.epsilon = float(L[1]) * u.pc
				
			elif L[0].lower() == "theta":
				self.theta = float(L[1])
				
			elif L[0].lower() == "time_bet_restartfiles":
				self.resfile = float(L[1])
				
			elif L[0].lower() == "time_bet_snapshots":
				self.snapfile = float(L[1]) * u.yr
		f.close()
				
	def check(self):
		try:
		
			A = self.epsilon
			A = self.theta
			A = self.t_end
			A = self.h
		except:
			raise ValueError("Error, important parameter missing from param file")
			
		try:
			A = self.resfile
		except:
			#set default
			self.resfile = 3600 * u.s
		try:
			A = self.snapfile
		except:
			self.snapfile = 1e7 * u.yr
				
class Node:
	def __init__(self , x , y , z , slength):
		self.x = x
		self.y = y
		self.z = z
		self.L = slength
		self.children = []
		self.particles = []
		self.com = None
	
	def in_node(self , particle):
		
		px = particle.r[t_ind][0]
		py = particle.r[t_ind][1]
		pz = particle.r[t_ind][2]
		
		if px < self.x + self.L and px >= self.x:
			if py < self.y + self.L and py >= self.y:
				if pz < self.z + self.L and pz >= self.z:
					return True
					
		return False
	

	def calc_all_com(self):
		'''
		computes the com for every node in our tree
		taakes in a time step t
		'''
		if len(self.particles) == 0:
			self.com = Vector( [ 0 * u.pc , 0 * u.pc , 0 * u.pc] )
			self.TM = 0 * u.Msun
			return 0 ###Just need to exit the function here
			
		elif len(self.particles) == 1:
			self.com = self.particles[0].r[t_ind]
			self.TM = self.particles[0].M
			return 0
		
		com = Vector( [ 0 * u.pc * u.M_sun , 0 * u.pc * u.M_sun, 0 * u.pc* u.M_sun] )
		Mass = 0 * u.M_sun
		
		for i in self.children:
			if i.com == None:
				i.calc_all_com()
			try:
				com += i.com * i.TM
			except:
				print (com.elements , i.com * i.TM)
				sys.exit()
			Mass += i.TM
		
		self.com = com * (1 / Mass)
		self.TM = Mass
			
		
	def reproduce(self):
	
		NL = self.L / 2.0
		LLL = Node(self.x , self.y , self.z , NL)
		LLU = Node(self.x , self.y , self.z + NL , NL)
		LUL = Node(self.x , self.y + NL , self.z , NL)
		ULL = Node(self.x + NL , self.y , self.z , NL)
		LUU = Node(self.x , self.y + NL , self.z + NL , NL)
		ULU = Node(self.x + NL , self.y , self.z + NL , NL)
		UUL = Node(self.x + NL , self.y + NL , self.z , NL)
		UUU = Node(self.x + NL , self.y + NL , self.z + NL , NL)
		###Now we need to place the existing particles into one of the child nodes.
		
		self.children = [LLL , LLU , LUL , ULL , LUU , ULU , UUL , UUU]
		for i in self.children:
			if i.in_node(self.particles[0]):
				i.add_particle(self.particles[0])
		
		
	def add_particle(self , particle):
		if len(self.particles) >= 1:
			self.particles.append(particle)
			if len(self.children) == 0:
				self.reproduce()
			
			
			for i in self.children:
			
				if i.in_node(particle):
				
					i.add_particle(particle)
					break
					
		
		elif len(self.particles) == 0:
		
			self.particles.append(particle)
		

class Vector:		

	def __init__(self , L):
		self.elements = L
		self.m = None
		
	def __getitem__(self , index):
		return self.elements[index]
		
	def __setitem__(self , index , val):
		self.elements[index] = val
	
	def __len__(self):
		return len(self.elements)
	
	def __add__(self , other):
	
		N = Vector([0] * len(self))
		for i in range(len(self)):
			N[i] = self[i] + other[i]
		return N
		
	def __sub__(self , other):
	
		N = Vector([0] * len(self))
		for i in range(len(self)):
			N[i] = self[i] - other[i]
			
		return N
		
	def __mul__(self , other):
		
		if type(other) == type(42) or type(other) == type(42.42) or type(other) == type(42.01 * u.pc):
			N = Vector( [0] * len(self))
			for i in range(len(N)):
				N[i] = self[i] * other
				
		elif type(other) == type(self):
			N = 0
			
			for i in range(len(self)):
				N += self[i] * other[i]
			
		return N
		
	def __eq__(self , other):
	
		if self.elements == other.elements:
			return True
		return False
		
		
	def transform(self , unit):
	
		for i in range(len(self.elements)):
			self.elements[i] = self.elements[i].to(unit)
	
	def print(self):
		print (self.elements)
		
		
	def unit(self):
		if self.m == None:
			self.mag()
		N = Vector([0] * len(self))
		for i in range(len(N)):
			N[i] = self[i] / self.m
		self.uv = N
		return N
		
		
	def mag(self):
		M = 0
		for i in self.elements:
			M += i ** 2
		self.m = np.sqrt(M)
		return self.m
		
	
class Particle:
	'''
	Class used to store information for a single particle
	takes in the x , y , z positions. Defaults to units of pc, unless inputs include astropy units
	Also takes in a mass, assumed to be in solar masses
	Can optionally take in a particle id
	'''
	
	def __init__(self , x , y , z , M , id):
		ct = type(5 * u.m)
		warn = False ##Triggers a unit warning if set to True
		if type(x) != ct:
		
			x *= u.pc
			warn = True
			
			
		if type(y) != ct:
			y *= u.pc
			warn = True
			
			
		if type(z) != ct:
			z *= u.pc
			warn = True
		if type(M) != ct:
			M *= u.Msun
			warn = True
			
			
		if warn:
			print ("Warning: At least on quantity lacked units, so the default was assumed")
		self.r =  [ Vector([x , y , z])]
		self.M = M
		self.id = id
		self.v = 0
		
		
	def __eq__(self , other):
		'''
		defines equality of particles. If id numbers are the same, then particles are the same
		if either particle lacks an id, then we use the starting position
		'''
		if self.id == other.id and self.id != None and other.id != None:
			return True
		return False
		
	def accel(self , other):
		"""
		Calculates the acceleration on this particle due to another particle
		Takes in self and another particle object
		Does not include force softening
		returns an acceleration vector
		"""
		
		
		
		rel = other.r[t_ind] - self.r[t_ind]
		
		###Magnitude of r12
		rel.mag()
		
		###r12^ (unit vector)
		rel.unit()
		a_mag = other.M * const.G / rel.m ** 2
		
		acc = rel.uv * a_mag
		return acc
		
	def soft_acc(self , other):
		
		"""
		Calculates the acceleration on this particle due to another particle
		Takes in self and another particle object
		Includes force softening
		returns an acceleration vector
		"""
	
		rel = other.r[t_ind] - self.r[t_ind]

		
		rel.mag()
		rel.unit()
		
		a_mag = other.M * const.G * S(rel) / rel.m
		
		acc = rel.uv * a_mag
		
		return acc
		
	def update_r(self , acc , h):
		'''
		Takes in an acceleration vector acc, step size h, and time step t_step
		updates the position of the particle based on this acceleration
		'''
		
		
		##Updates r given an acceleration and step size and t_step
		nr = copy.deepcopy(self.r[t_ind])
		nr *= 2
		nr -= self.r[0]
		#nr += nr + acc * h ** 2
		nr += acc * h ** 2
		self.r.append(nr)
		
	def ndt_update_r(self, acc , old_h , new_h):
	
		##Updates r given an acceleration and step size and t_step
		nr = copy.deepcopy(self.r[t_ind])
		nr += (nr - self.r[0]) * (new_h / old_h)
		nr += acc * new_h ** 2

		self.r.append(nr)
		
	def distance(self , other):
		#Complutes the distance between two particles at a time t
		r = self.r[t_ind] - other.r[t_ind]
		return r.mag()

class Sim:
	def __init__(self , Parts):
		self.particles = Parts
	def __getitem__(self , index):
		return self.particles[index]
	def __eq__(self , other):
		if self.particles == other.particles:
			return True
		return False
		
	def write_snapshot(self , fname):
		np.save(fname , self)
		
	def read_snapshot(self , fname):
		self.particles = np.load(fname , allow_pickle = True)[()].particles
		
	def write_restart_files(self , ctime):
		##TODO finish writing restart file code
		try:
			shutil.rmtree("restartfiles")
		except:
			print ("Creating first restart file")
		os.mkdir("restartfiles")
		self.write_snapshot("restartfiles/restart")
		shutil.copyfile(param.paramfile , "restartfiles/params.txt")
		f2 = [ctime.to(u.yr)  , param.paramfile]
		np.save("restartfiles/par" , f2)
		
		
	def clear(self):
		for i in range(len(self.particles)):\
			del self.particles[i].r[0]
	
	def make_plot(self):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		X  = []
		Y = []
		Z = []
		for i in self.particles:
			X.append(i.r[1][0].to("Mpc").value)
			Y.append(i.r[1][1].to("Mpc").value)
			Z.append(i.r[1][2].to("Mpc").value)
		ax.set_xlabel('X (Mpc)')
		ax.set_ylabel('Y (Mpc)')
		ax.set_zlabel('Z (Mpc)')
		ax.scatter(X, Y, Z, zdir='z')	
		plt.show()
	
			

def direct_summation(t_end , h , IC):
	'''
	This is an nbody solver that uses direct summation
	Not useful for our problem, but it is good for testing purposes
	takes in an end time and step size, t_end and h
	takes in a list of particles, IC
	returns a new list of particles
	'''
	
	t = 0 * t_end ##Preserves units
	t_step = 1
	while t < t_end:
		for P1 in IC:
			total_acc = Vector([0 , 0 , 0])
			for P2 in IC:
				if P1 == P2:
					continue
				acc = P1.accel(P2)
				total_acc += acc
			P1.update_r(total_acc , h)
		t_step += 1
		t += h
	return IC
	
	
def Find_Tree(IC):
	'''
	This takes in a list of particle objects, IC, and a time t.
	t should be an integer describing which time step we are interested in
	This function generates and returns a tree
	returns a node object
	'''
	x = []
	
	y = []
	
	z = []
	for i in IC:
	

		x.append(i.r[t_ind][0])
		y.append(i.r[t_ind][1])
		z.append(i.r[t_ind][2])
	
	lunit = i.r[t_ind][0].unit
	
	sx = min(x) - 1 * lunit
	sy = min(y) - 1 * lunit
	sz = min(z) - 1 * lunit
	
	L = max([ max(x) - sx + 1 * lunit , max(y) - sy+ 1 * lunit , max(z) - sz + 1 * lunit])
	Pnode = Node(sx , sy , sz , L)
	
	
	for i in IC:
	
		Pnode.add_particle( i)
	
	
	Pnode.calc_all_com()
	return Pnode
	
def BH_Acceleration(IC  , Tree , Part):

	'''
	Calculates an acceleration using the Barnes - Hut tree algorithm
	Takes in a list of particles IC and time step t
	Takes in a tree (Node object) for the particle list
	Part is a particle object
	returns the acceleration vector for Part
	'''
	
	
	Acc = Vector ( [ 0 , 0 , 0 ] )
	L = Tree.L
	TP = Particle(Tree.com[0] , Tree.com[1] , Tree.com[2] , 1 * u.Msun , -1)
	TP.r.append(Vector([Tree.com[0] , Tree.com[1] , Tree.com[2]]))
	D = Part.distance(TP)
	if D.value == 0:
		return Acc
	T = L / D
	if len(Tree.children) == 0 and len(Tree.particles) == 1:
		if Tree.particles[0] == Part:
			return Acc
			
		else:
			
			Acc = Part.soft_acc(Tree.particles[0])
			
	elif T < param.theta:
	
		
		Tree_part = Particle(  Tree.com[0] , Tree.com[1] , Tree.com[2] , Tree.TM , -1)
		Tree_part.r.append(Vector([ Tree.com[0] , Tree.com[1] , Tree.com[2]]))
		Acc += Part.soft_acc(Tree_part)
		
	elif len(Tree.children) > 0:
		for i in Tree.children:
			Acc += BH_Acceleration(IC , i , Part)
	return Acc
	
def Barnes_Hut(IC , ctime = 0):

	h = param.h
	t_end = param.t_end
	'''
	Given a list of particles, will run an n-body simulation using the Barnes - Hut tree algorith
	IC is a list of particles
	will return the list of particles complete with updated positions
	
	'''
	
	t = ctime * u.yr ##Preserves units
	out = 1
	ts = 1
	start = time.time()
	tbs = int(param.snapfile / h)
	if tbs == 0:
		tbs = q
	while t < t_end:
		Tree = Find_Tree(IC)
		
		for P1 in IC:
			a = BH_Acceleration(IC , Tree , P1)
			P1.update_r(a , h)
		IC.clear()
		
		if ts == tbs:
			IC.write_snapshot("snapshot_" + str(int(t)))
			out += 1
			ts = 1
			
		if (time.time() - start) >= param.resfile * 60:
			IC.write_restart_files(t)
			start = time.time()
			
		t += h
	return IC
	
def restart():
	global param
	param = Params("restartfiles/params.txt")
	A = np.load("restartfiles/restart.npy" , allow_pickle = True)
	IC = Sim(A)
	B = np.load("restartfiles/par.npy" , allow_pickle = True)
	Res = Barnes_Hut(IC , ctime = B[0])
	return Res
	

def change_dt(IC , h , nh , t_end , t_step = 1):
	'''
	This function will run a simulation given two snapshots with delta_t = h
	Will run a simulation with timestep nh
	'''
	
	###Step 1
	if h == nh:
		Res = Barnes_Hut(IC , nh , t_end)
		return Res
	else:
		##Take one step to get new time step
		Tree = Find_Tree(IC)
		
		for P1 in IC:
			a = BH_Acceleration(IC , Tree , P1)
			P1.ndt_update_r(a , h , nh)
		
		for i in range(len(IC)):
			IC[i].r.pop(1)
		
		Res = Barnes_Hut(IC , nh)
		return Res
		
		
	

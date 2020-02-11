import astropy.units as u
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import astropy.constants as const
import calc
import matrix

#rc('text', usetex=True)
plt.rc('font', size=14)

def vc(r , r200 , v200 , c):
	'''
	takes in the radius, r200 , v200 , and concentration for a halo
	returns the circular velocity at the radius r
	'''
	x = r / (v200.value * u.kpc)
	a = ((np.log(1 + c * x) - (c * x) / (1 + c * x)) / (np.log(1 + c) - c / (1 + c)))
	return np.sqrt((a / x) * v200 ** 2)
	
def Menc(r , r200 , v200 , c):

	'''
	takes in the radius, r200 , v200 , and concentration for a halo
	returns the total mass enclosed within the radius r
	'''
	
	v_c = vc(r  , r200 , v200 , c)
	return (r * v_c ** 2 / const.G).to(u.solMass)
	
def rho(r , r200 , v200 , c):
	'''
	takes in the radius, r200 , v200 , and concentration for a halo
	returns the density at a radius r
	'''
	
	rs = r200 / c
	return 1 / ((r / rs) * ((1 + r / rs) ** 2))
	
def M(r , r200 , v200 , c):
	'''
	takes in the radius, r200 , v200 , and concentration for a halo
	returns the mass inside of a spherical shell
	'''
	
	return 4 * np.pi * rho(r , r200 , v200 , c) * r ** 2

def J(T , l , up):

	'''
	Takes in a temperature, and two energy states (which both should be integers)
	returns the mean intensity of the radiation field at the frequency corresponding tho the energy of this transition
	Assumes we are dealing with hydrogen
	'''
	
	
	nu = get_freq(l , up).to(1 / u.s)
	C = (2 * const.h * ((nu) ** 3)) / (const.c ** 2)
	return C / (np.exp((const.h * nu / (const.k_B * T)).to(u.dimensionless_unscaled)) - 1)
	
def get_A(l , up):

	'''
	reads in our table of Einstein coefficients
	returns the value of A for a given set of energy levels l and up
	'''
	
	
	in_file = open("A_coefficients.dat")
	A = matrix.Matrix((9 , 9))
	for i in in_file.readlines():
		if i[0] == "#":
			continue
		a = int(i.split()[0].strip().replace("," , ""))
		b = int(i.split()[1].strip().replace("," , ""))
		A.elements[a - 1][b - 1] = float(i.split()[2].strip().replace("," , ""))
	in_file.close()
	return A.get(l , up) / u.s
	
def get_B(l , up):
	'''
	calculates the Einstein coefficient called B
	'''
	Alu = get_A(l , up)
	a = (const.c ** 2) / (2 * const.h * get_freq(l , up).to(1 / u.s) ** 3)
	return a * Alu


def get_freq(l , up):
	'''
	finds the frequency for a particular line transition of hydrogen
	takes in two energy levels (integers)
	returns a frequency
	'''
	E = abs(( -13.6 * u.eV / (up ** 2) ) - ( -13.6 * u.eV / (l ** 2) ))
	return E / const.h
	
def densities(T):

	'''
	This function takes in a temperature T
	computes all of our number densities
	returns a Matrix object containing all of the number densities
	'''
	
	###(3x3)
	b = matrix.Matrix((3 , 1) , ([[0] , [0], [1]]))
	A = matrix.Matrix((3 , 3))
	a_00 = ((get_B(1 , 2) * J(T , 1 , 2) + get_B(1 , 3) * J(T , 1 , 3))).value
	a_02 = -1 * get_A(1 , 3) + get_B(1 , 3) * J(T , 1 , 3)
	A.elements[0][0] = a_00
	A.elements[0][1] = -1 * (get_A(1 , 2) + get_B(1 , 2) * J(T , 1 , 2)).value
	A.elements[0][2] = -1 * (get_A(1 , 3) + get_B(1 , 3) * J(T , 1 , 3)).value
	A.elements[1][1] = (get_B(2 , 1) * J(T , 2 , 1) + get_A(2 , 1) + get_B(2 , 3) * J(T , 2 , 3)).value
	A.elements[1][0] = -1 * (get_B(1 , 2) * J(T  , 1 , 2)).value
	A.elements[1][2] = -1 * (get_B(2 , 3) * J(T , 2  ,3) + get_A(2 , 3)).value
	'''
	A.elements[2][0] = -1 * (get_B(1 , 3) * J(T , 1 , 3)).value
	A.elements[2][1] = -1 * (get_B(2 , 3) * J(T , 2 , 3)).value
	A.elements[2][2] = (get_B(3 , 1) * J(T , 3 , 1) + get_A(3 , 1) + get_A(3, 2) + get_B(3 , 2) * J(T , 3 , 2)).value
	'''
	A.elements[2][0] = 1
	A.elements[2][1] = 1
	A.elements[2][2] = 1
	x = matrix.solve_eq(A , b)
	return x


###Runs code to produce all of our plots

r200 = 200 * u.kpc
v200 = 200 * u.km / (u.s)
c = 9.39
x = range(1 , 500)
y = []
yp = []
mint = []
f = calc.centered_diff(lambda x: M(x * u.kpc , r200 , v200 , c).value , 1e-5)
for i in x:
	y.append(np.log10(M(i * u.kpc , r200 , v200 , c).value))
	yp.append(f(i))
	mint.append(np.log10(Menc(i * u.kpc , r200 , v200 , c).value))
plt.plot(x , y)
plt.xlabel("r (kpc)")
plt.ylabel("M(r) M_{sun}")
plt.show()

plt.plot(x , yp)
plt.xlabel("r (kpc)")
plt.ylabel("M'(r)")
plt.show()


plt.plot(x , mint)
plt.xlabel("r (kpc)")
plt.ylabel("M_enc")
plt.show()

	
Mtot = Menc(1e8 * u.kpc , r200 , v200 , c)
print (np.log10(Mtot.value))

print (M(5 * u.kpc , r200 , v200 , c))
	
	
###On to problem 5
N = 1 / (u.cm ** 3)
T_arr = list(range(1000 , 20000 , 100))
y = []
'''
for i in T_arr:
	y.append(densities(i * u.K).elements[0][0])
print (T_arr , y)
plt.plot(T_arr , y)
plt.xlabel("T (K)")
plt.ylabel("Density")
plt.show()

'''

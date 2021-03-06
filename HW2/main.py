import astropy.units as u
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import astropy.constants as const
import calc
import matrix
import copy
rc('text', usetex=True)
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
	
def get_A(upper , lower):

	'''
	reads in our table of Einstein coefficients
	returns the value of A for a given set of energy levels l and up
	'''
	
	
	in_file = open("A_coefficients.dat")
	A = matrix.Matrix((9 , 9))
	for i in in_file.readlines():
		if i[0] == "#":
			continue
		up = int(i.split()[0].strip().replace("," , ""))
		low = int(i.split()[1].strip().replace("," , ""))
		A.elements[up - 1][low - 1] = float(i.split()[2].strip().replace("," , ""))
	in_file.close()
	return A.get(lower , upper) / u.s
	
def get_B(upper , lower):
	'''
	calculates the Einstein coefficient called B
	'''
	Alu = get_A(upper , lower)
	if Alu == 0:
		Alu = (get_A(lower , upper)) * ( ((2 * lower ** 2) / (2 * upper ** 2)) ** 1)
		
	
	a = (const.c ** 2) / (2 * const.h * get_freq(lower , upper).to(1 / u.s) ** 3)
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
	
	A.elements[0][0] = ((get_B(1 , 2) * J(T , 1 , 2) + get_B(1 , 3) * J(T , 1 , 3))).value
	A.elements[0][1] = -1 * (get_A(2 , 1) + get_B(2 , 1) * J(T , 2 , 1)).value
	A.elements[0][2] = -1 * (get_A(3 , 1) + get_B(3 , 1) * J(T , 3 , 1)).value
	A.elements[1][1] = (get_B(2 , 1) * J(T , 2 , 1) + get_A(2 , 1) + get_B(2 , 3) * J(T , 2 , 3)).value
	A.elements[1][0] = -1 * (get_B(1 , 2) * J(T  , 1 , 2)).value
	A.elements[1][2] = -1 * (get_B(3 , 2) * J(T , 3  ,2) + get_A(3 , 2)).value
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


def nden(T):
	'''
	This function takes in a temperature T
	computes all of our number densities
	returns a Matrix object containing all of the number densities
	'''
	
	n_levels = 9
	
	A = matrix.Matrix((n_levels , n_levels))
	b = matrix.Matrix((n_levels , 1))
	b.elements[n_levels - 1][0] = 1
	
	for i in range(n_levels):
		###These for loops set up the matrix A, which contains many einstein coefficieints
		### A contains all of the coefficients for a system of equations for the various number densities
		
		for j in range(n_levels):
			if i == n_levels - 1:
				A.elements[i][j] = 1
				continue
			elm = 0
			if i == j:
				for k in range(n_levels):
					if k == i: ###No transitions from one state to itself
						continue
					elm += (get_A(i + 1 , k + 1).value)
					elm += (get_B(i + 1 , k + 1) * J(T , i + 1 , k + 1)).value
				A.elements[i][j] = elm
				
			else:
				#for k in range(n_levels):

				elm = (get_B(j + 1 , i + 1) * J(T , j + 1, i + 1) + get_A(j + 1, i + 1)).value
				A.elements[i][j] = -1 * elm
		

	###Now we solve our matrix problem
	
	x = matrix.solve_eq(A , b)

	return x
###Runs code to produce all of our plots


##Sets a few parameters
r200 = 200 * u.kpc
v200 = 200 * u.km / (u.s)
c = 9.39

####

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
plt.ylabel("$M (M_{\odot})$")
plt.savefig("P21.pdf")
#plt.show()
plt.close()

plt.plot(x , yp)
plt.xlabel("r (kpc)")
plt.ylabel("$dM(r)/dr (M_{\odot})$")
plt.savefig("P22.pdf")
plt.close()

plt.plot(x , mint)
plt.xlabel("r (kpc)")
plt.ylabel("$M_{enc} (M_{\odot})$")
plt.savefig("P23.pdf")
plt.close()

	
Mtot = Menc(1e8 * u.kpc , r200 , v200 , c)
print (np.log10(Mtot.value))

print (M(5 * u.kpc , r200 , v200 , c))
	
	
###On to problem 5
N = 1 / (u.cm ** 3)


T = 1e4
y1 = []
y2 = []
y3 = []
y4 = []
y5 = []
y6 = []
y7 = []
y8 = []
y9 = []

T_arr = []
while T < 1e5:
	T_arr.append(np.log10(T))
	X = nden(T * u.K)
	y1.append(X.elements[0][0])
	y2.append(X.elements[1][0])
	y3.append(X.elements[2][0])
	y4.append(X.elements[3][0])
	y5.append(X.elements[4][0])
	y6.append(X.elements[5][0])
	y7.append(X.elements[6][0])
	y8.append(X.elements[7][0])
	y9.append(X.elements[8][0])
	T *= 1.1
	
f , ax = plt.subplots(3 , 3 , sharex = True , sharey = True , figsize = (10 , 10))

plot_u = "$( cm^{-3} )$"
ax[0,0].plot(T_arr , y1)
ax[0,0].set_xlabel("$log(T) (K)$")
ax[0,0].set_ylabel("$n_1$" + plot_u)


ax[0,1].plot(T_arr , y2)
ax[0,1].set_xlabel("$log(T) (K)$")
ax[0,1].set_ylabel("$n_2$" + plot_u)

ax[0,2].plot(T_arr , y3)
ax[0,2].set_xlabel("$log(T) (K)$")
ax[0,2].set_ylabel("$n_3$" + plot_u)


ax[1,0].plot(T_arr , y4)
ax[1,0].set_xlabel("$log(T) (K)$")
ax[1,0].set_ylabel("$n_4$" + plot_u)

ax[1,1].plot(T_arr , y5)
ax[1,1].set_xlabel("$log(T) (K)$")
ax[1,1].set_ylabel("$n_5$" + plot_u)

ax[1,2].plot(T_arr , y6)
ax[1,2].set_xlabel("$log(T) (K)$")
ax[1,2].set_ylabel("$n_6$" + plot_u)

ax[2,0].plot(T_arr , y7)
ax[2,0].set_xlabel("$log(T) (K)$")
ax[2,0].set_ylabel("$n_7$" + plot_u)

ax[2,1].plot(T_arr , y8)
ax[2,1].set_xlabel("$log(T) (K)$")
ax[2,1].set_ylabel("$n_8$" + plot_u)

ax[2,2].plot(T_arr , y9)
ax[2,2].set_xlabel("$log(T) (K)$")
ax[2,2].set_ylabel("$n_9$"+ plot_u)

plt.subplots_adjust(hspace = .31)
plt.savefig("P5.pdf")
plt.close()

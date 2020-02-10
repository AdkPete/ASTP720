import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import calc
import matrix


def vc(r , v200 , c):
	x = r / (v200.value * u.kpc)
	a = ((np.log(1 + c * x) - (c * x) / (1 + c * x)) / (np.log(1 + c) - c / (1 + c)))
	return np.sqrt((a / x) * v200 ** 2)
	
def Menc(r , v200 , c):
	v_c = vc(r  ,v200 , c)
	return (r * v_c ** 2 / const.G).to(u.solMass)
	
def J(T):
	
	2 * const.h
	return T
	
in_file = open("A_coefficients.dat")

	
v200 = 200 * u.km / (u.s)
c = 9.39
x = range(1 , 300)
y = []
for i in x:
	y.append(np.log10(Menc(i * u.kpc , v200 , c).value))
	
plt.plot(x , y)
#plt.show()
	
M = Menc(1e8 * u.kpc , v200 , c)
print (np.log10(M.value))
	
N = 1 / (u.cm ** 3)

A = matrix.Matrix((9 , 9))
for i in in_file.readlines():
	if i[0] == "#":
		continue
	l = int(i.split()[0].strip().replace("," , ""))
	u = int(i.split()[1].strip().replace("," , ""))
	A.elements[l - 1][u - 1] = float(i.split()[2].strip().replace("," , ""))
	
	
A_12 = A.get(1 , 2)

M = matrix.Matrix((3 , 3))


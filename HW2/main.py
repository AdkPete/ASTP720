import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import calc
import matrix


def vc(r , r200 , v200 , c):
	x = r / (v200.value * u.kpc)
	a = ((np.log(1 + c * x) - (c * x) / (1 + c * x)) / (np.log(1 + c) - c / (1 + c)))
	return np.sqrt((a / x) * v200 ** 2)
	
def Menc(r , r200 , v200 , c):
	v_c = vc(r  , r200 , v200 , c)
	return (r * v_c ** 2 / const.G).to(u.solMass)
	
def rho(r , r200 , v200 , c):
	
	rs = r200 / c
	return 1 / ((r / rs) * ((1 + r / rs) ** 2))
	
def M(r , r200 , v200 , c):
	return 4 * np.pi * rho(r , r200 , v200 , c) * r ** 2

def J(T , l , up):
	nu = get_freq(l , up)
	C = 2 * const.h *(nu) ** 3 / const.c ** 2
	print ((const.h * nu / (const.k_B * T)).to(u.dimensionless_unscaled))
	return C / (np.exp((const.h * nu / (const.k_B * T)).to(u.dimensionless_unscaled)) - 1)
	
def get_A(l , up):
	A = matrix.Matrix((9 , 9))
	for i in in_file.readlines():
		if i[0] == "#":
			continue
		l = int(i.split()[0].strip().replace("," , ""))
		u = int(i.split()[1].strip().replace("," , ""))
		A.elements[l - 1][up - 1] = float(i.split()[2].strip().replace("," , ""))
		
	return A.get(l , up)
	
def get_B(l , up):
	Alu = get_A(l , up)
	return const.c ** 2 * Alu / (2 * const.h * get_freq(l , up))


def get_freq(l , up):
	E = abs(( -13.6 * u.eV / (up ** 2) ) - ( -13.6 * u.eV / (l ** 2) ))
	return E / const.h
	
in_file = open("A_coefficients.dat")

r200 = 230 * u.kpc
v200 = 200 * u.km / (u.s)
c = 9.39
x = range(1 , 300)
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
plt.ylabel("M(r) M_sun")
plt.show()

plt.plot(x , yp)
plt.xlabel("r (kpc)")
plt.ylabel("M'(r) M_sun")
plt.show()

plt.plot(x , mint)
plt.xlabel("r (kpc)")
plt.ylabel("M_enc (r) M_sun")
plt.show()

	
Mtot = Menc(1e8 * u.kpc , r200 , v200 , c)
print (np.log10(Mtot.value))

print (M(5 * u.kpc , r200 , v200 , c))
	
	
###On to problem 5
N = 1 / (u.cm ** 3)




import numpy as np
import astropy.units as u
from galpy.potential import MWPotential2014
from galpy.potential import evaluateRforces
from galpy.potential import evaluatezforces
from galpy.util import bovy_conversion
import matplotlib.pyplot as plt
from genetic_algorithms import Genetic
import scipy.stats as stats

poten = MWPotential2014

class Vector:

	def __init__(self , x , y , z):
		
		self.x = x
		self.y = y
		self.z = z
		
		self.mag = np.sqrt(x ** 2 + y ** 2 + z ** 2)
		
	def dot(self , other):
		return self.x * other.x + self.y * other.y + self.z * other.z
	
	def scale(self , val):
		self.x *= val
		self.y *= val
		self.z *= val
		
	
	def project(self , other):
		##Will project other onto self
		proj = Vector(self.x , self.y , self.z)
		D = self.dot(other)
		D /= (self.mag ** 2)
		proj.scale(D)
		
		return Vector(proj.x , proj.y , proj.z)

def rv(x , y , z):
	
	conv = bovy_conversion.force_in_pcMyr2(220.,8) *  u.pc / (u.Myr ** 2) ##Unit conversion
	
	x *= u.kpc
	y *= u.kpc
	z *= u.kpc
	
	##Sun's location
	sx = 8 * u.kpc 
	sy = 0 * u.kpc
	sz = 4 * u.kpc
	
	sun_ar = evaluateRforces(poten , sx , sz , phi = 0 * u.deg) * conv
	sun_az = evaluatezforces(poten , sx , sz , phi = 0 * u.deg) * conv
	
	
	sun_acc = Vector(sun_ar , 0 * conv , sun_az)
	
	star_ar = evaluateRforces(poten , np.sqrt(x ** 2 + y ** 2) , z , phi = 0 * u.deg ) * conv
	star_az = evaluatezforces(poten , np.sqrt(x ** 2 + y ** 2) , z , phi = 0 * u.deg ) * conv
	
	star_acc = Vector(star_ar , 0 * conv , star_az)
	
	
	heliocentric_acc = Vector(star_acc.x - sun_acc.x , star_acc.y - sun_acc.y , star_acc.z - sun_acc.z)
	
	r = Vector(sx - x , sy - y , sz - z )
	
	acc = r.project(heliocentric_acc)
	
	

	delta_v = (acc.mag * 10 * u.yr).to(u.cm / u.s)
	
	return delta_v , r
	
def f(X):
	
	dv , r = rv(X[0] , X[1] , X[2])
	
	factor = stats.norm.pdf(r.mag.value , 0 , 1.5)
	
	return dv.value * factor
	
opt = Genetic(f , [15 , 15] , [ [ 3 , 20 ] , [-5 , 5 ] , [-10 , 10 ] ] , 1000)
opt.initialize()

opt.update(250)

print (opt.creatures[0].get_params() , opt.creatures[0].fitness)

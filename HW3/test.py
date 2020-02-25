import ode

import numpy as np
import unittest

def example_func(t , y):
	return t ** 2 + 5 * t + 3
	
def example_sltn(t):
	return (1.0 / 3) * t ** 3 + (5.0 / 2) * t ** 2 + 3 * t + 7
	
def example_solve():
	
	##solves a simple ode using all three methods
	##this is an example of how to use solve_ode
	
	E1 = ode.solve_ode(example_func , 0 , 7 , t_end = 5) ###Sets up a solver with the default step size h = 1e-3 from 0 to 5
	
	###Solve with all three ode solvers
	
	E1.Forward_Euler()
	
	print (E1.fe_t[-1] , E1.fe_y[-1])
	
	E1.Heun()
	
	print (E1.heun_t[-1] , E1.heun_y[-1])
	
	E1.RK4()
	
	print (E1.rk4_t[-1] , E1.rk4_y[-1])
	
	###We can also solve this problem specifying a number of steps instead of a step size
	
	E2 = ode.solve_ode(example_func , 0 , 7 , t_end = 5 , n = 1e4) ##Uses 10000 steps across the interval
	
	E2.Forward_Euler()
	
	E2.Heun()
	
	E2.RK4()
	
class TestODE(unittest.TestCase):

	###Tests that our three methods give reasonable results. Note that it is required that Forward Euler has 10 times as many steps to reach the same accuracy
	
	def test_fe(self):
		Test = ode.solve_ode(example_func , 0 , 7 , t_end = 5 , n = 2e5)
		Test.Forward_Euler()
		
		self.assertAlmostEqual(round(Test.fe_y[-1] , 3) , round(example_sltn(5) , 3))
	
	def test_heun(self):
		Test = ode.solve_ode(example_func , 0 , 7 , t_end = 5 , n = 2e4)
		Test.Heun()
		
		self.assertAlmostEqual(round(Test.heun_y[-1] , 3) , round(example_sltn(5) , 3))
		
		
	def test_rk4(self):
		
		Test = ode.solve_ode(example_func , 0 , 7 , t_end = 5 , n = 2e4)
		Test.RK4()
		
		self.assertAlmostEqual(round(Test.rk4_y[-1] , 3) , round(example_sltn(5) , 3))
		

#example_solve()
unittest.main()

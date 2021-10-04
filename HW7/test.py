import unittest
import time
import matplotlib.pyplot as plt
from FFT import *

class Test_All(unittest.TestCase):
	
	def test_slow_fft(self):
		x = np.linspace(0 , 10 , 100)
		t = np.linspace(0 , 10 , 100)
		x = np.sin(2 * np.pi * x)
		X = slow_dft(x)
		
		NX = np.fft.fft(x)
		
		self.assertEqual(np.allclose(X , NX) , True)
		
	def test_fast_fft(self):
		
		x = np.linspace(0 , 10 , 256)
		t = np.linspace(0 , 10 , 256)
		x = np.sin(2 * np.pi * x)
		X = fft(x)
		
		
		NX = np.fft.fft(x)
		
		self.assertEqual(np.allclose(X , NX) , True)

		
	def test_speedup(self):
		x = np.linspace(0 , 10 , 1024)
		t = np.linspace(0 , 10 , 1024)
		x = np.sin(2 * np.pi * x)
		
		start = time.time()
		X = slow_dft(x)
		slow_time = time.time() - start
		
		
		start = time.time()
		X = fft(x)
		fast_time = time.time() - start
		
		self.assertGreater(slow_time/ 10 , fast_time)
		
	def test_fast_ifft(self):
		
		x = np.linspace(0 , 10 , 256)
		t = np.linspace(0 , 10 , 256)
		x = np.sin(2 * np.pi * x)
		X = fft(x)
		
		NX = ifft(X)
		NPX = np.fft.ifft(X)

		
		self.assertEqual(np.allclose(NPX , NX) , True)
		
if __name__ == "__main__":
	unittest.main()

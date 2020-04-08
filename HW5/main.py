from least_squares import *
import numpy as np
def read_file(fname):
	f = open(fname)
	Z = []
	M = []
	P = []
	for i in f.readlines():
		
		if i[0] == "#":
			continue
		else:
			P.append(np.log10(float(i.split(",")[1])))
			m = float(i.split(",")[3])
			distance = float(i.split(",")[2])
			d = distance = 1000.0 ##pc
			E = float(i.split(",")[-2])
			Rv = 3.1
			Av = Rv * E
			
			AM = m - 5 * np.log10(d) + 5 + Av
			M.append(AM)
			Z.append(float(i.split(",")[-1]))
	return P , M , Z
		
LP , M , Z = read_file("cepheid_data.txt")
alpha , beta , gamma = more_general_fit(M , LP , Z)
mk_3d(M , LP , Z , alpha , beta , gamma)
print (alpha , beta , gamma)

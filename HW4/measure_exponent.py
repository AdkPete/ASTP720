import numpy as np
import nbody
import sys

F1 = sys.argv[1]
F2 = sys.argv[2]


R1 = nbody.read_snapshot(F1)
R2 = nbody.read_snapshot(F2)
tid = 908
for i in R1:
	if i.id == tid:
		i1r = i.r[0]
		f1r = i.r[-1]
		

tid = 908
for i in R2:
	if i.id == tid:
		i2r = i.r[0]
		f2r = i.r[-1]

print ((i1r - i2r).mag())
print ((f1r - f2r).mag())

d1 = (i1r - i2r).mag()
d2 = (f1r - f2r).mag()


### d2 = e^(lambda t) d1

Lt = np.log((d2 / d1).value)
L = Lt / 1
print (L)

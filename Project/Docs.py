
## THis script has some examples / nots about how to use my code


## The example here is a simple one, finding the max value of a gaussian pdf
## In thise case we will use a 4-D PDF, with a peak value at (0 , 0 , 0 , 0)
## It also has varied standard deviations, although this algorithm will not care
## The lack of depencece on derivatives is a nice feature of these algorithms
## This can lead to excellent performance relative to other methods on some problems


import genetic_algorithms as gene
import scipy.stats as stats

def test_f(X):

	## 3d fitness function
	sig = [ [1.5 , 0 , 0 , 0], [0 , 74 , 0 , 0 ] , [0 , 0 , 100 , 0] , [0 , 0 , 0 , 7] ]
	mu = [ 0 , 0  , 0 , 0]

	rv = stats.multivariate_normal( mu , sig )

	return rv.pdf(X)	
	
## So the function we want to maximize is test_f
## Now let us assume that we know almost nothing about the right answer
## As such, we set the following area to search in, from -500 to 500 in all directions

B = [-500 , 500]

bound = [B , B , B , B]

##Now we can set up the algorithm

opt = gene.Genetic(test_f , bound) ##Defaults should work well here

## Run for 10 iterations, return best result

opt.update(10)

print (opt.creatures[0].get_params())

##This probably didn't quite get there, so we can continue for some more genetations

opt.update(50)

print (opt.creatures[0].get_params())

## That should get pretty close

## We could do better with a larger pop size, smaller bounds or more iterations


## If you need a more precise answer, increase the precision
## Higher precision is slower, but it determines how long a DNA string is
## Longer DNA lets you represent more numbers, but is slower to work with

new = gene.Genetic(test_f , bound , popsize = 5000 , precision = [35 , 35])

##Note that this is uninitialized, so contains zero creatures
## can be set up with new.initialize() and new.determine_fitness()
## Alternatively new.update() will do it automatically

## This would be slow because of the larger pop size and higher precision
## but should give a good answer


new.update(5)


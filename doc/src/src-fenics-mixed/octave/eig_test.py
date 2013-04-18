import numpy, pytave
A = numpy.random.random([2,2])
e = pytave.feval(1, "eig", A)
print 'Eigenvalues:', e

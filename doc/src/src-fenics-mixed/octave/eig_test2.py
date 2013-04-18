import numpy, pytave
A = numpy.random.random([2,2])
B = numpy.random.random([2,2])
e, v = pytave.feval(2, "eig", A, B)
print 'Eigenvalues:', e
print 'Eigenvectors:', v



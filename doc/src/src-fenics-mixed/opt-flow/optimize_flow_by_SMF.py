import numpy
import pytave
import sys

def dump_vector (filename, name, AA):
    """Dump vector to matlab readable format."""
    f = open(filename, 'w')
    for i in range(AA.shape[0]):
        if abs(AA[i]) > 10e-10:
            f.write("%s (%d) = %e;\n " % (name,i+1,AA[i]))

def dump_matrix(filename, name, AA):
    """Dump matrix to matlab readable format."""
    f = open(filename, 'w')
    for i in range(AA.shape[0]):
        for j in range(AA.shape[1]):
           if abs(AA[i,j]) > 10e-10:
               f.write("%s (%d, %d) = %e;\n " %
                       (name,i+1,j+1,AA[i,j]))

def monitor(Ai, J, nit, curr_bestA, curr_bestJ,
            delta, prev_it, improve, spc):
    """Monitor convergence in Ai*.m, J*.m, and output*.m."""
    dump_matrix("Ai_%d.m" % nit , "Ai_%d" % nit, Ai)
    dump_vector("J_%d.m" % nit , "J_%d" % nit, J)
    f = open("output_%d.tmp"% nit, "w")
    f.write("nit %d\n" % nit)
    f.write("curr_bestA %e, %e, %e\n" %
            (curr_bestA[0], curr_bestA[1], curr_bestA[2]) )
    f.write("curr_bestJ %e\n" % curr_bestJ)
    f.write("delta %e\n" % delta)
    f.write("prev_it %s\n" % prev_it)
    f.write("improve %d\n" % int(improve))
    f.write("spc %d, %d, %d\n" % (spc[0][0], spc[0][1], spc[0][2]))
    f.close()

def run_simulations(Ai):
    """Run a sequence of simulations with input parameters Ai."""
    import flow_problem
    plot = True
    if len(Ai.shape) == 1: # only one set of parameters
        J = numpy.zeros(1)
	K0, x0, c = Ai
	p = flow_problem.FlowProblem2Optimize(K0, x0, c, plot)
	goal1, goal2 = p.run()
	J[0] =  goal1 + goal2
    else: # several sets of parameters
        J = numpy.zeros(len(Ai))
        i = 0
	for a in Ai:
	    K0, x0, c = a
	    p = flow_problem.FlowProblem2Optimize(K0, x0, c, plot)
	    goal1, goal2 = p.run()
	    J[i] = goal1 + goal2
	    i = i+1
    return Ai, J

def initialize():
    """Create initial sampling points."""
    AiLHS = pytave.feval(1, "initial_pts_Nd", ninit, amin, amax, N)
    AiLHS = numpy.array(AiLHS)
    Ai = numpy.zeros([ninit,3])
    for i in range(0,ninit):
	x  = pytave.feval(1, "find_near_pt",
                          AiLHS[0, i,:], N, spc, amin, amax)
        x = numpy.array(x)
	Ai[i,:] = x[0][0]
    return Ai

def coarsen(Ai):
    # not implemented yet
    return Ai

def check_for_new_points(next_ptsall, Aall):
    """Check that the new points have not been found before."""
    next_pts = None
    for i in range(0,len(next_ptsall)):
        match = 0
        for j in range(Aall.shape[0]):
    	    d = next_ptsall[i,:]-Aall[j,:]
            if max(abs(d)) <=1e-10:
	        match = 1
        if match == 0:
             if isinstance(next_pts, numpy.ndarray):
	         next_pts = numpy.vstack((next_pts,next_ptsall[i,:]))
             else:
	         next_pts = next_ptsall[i,:]
    return next_pts


def search(Aall, Jall, curr_bestA, theta, upb,
           lob, N, amin, amax, spc, delta):
    """Search step."""
     # make sure that all points are unique
    (Am, Jm) = pytave.feval(2, "dsmerge", Aall, Jall)

    next_ptsall = []
    next_pts = None
    max_no_searches = 100
    no_searches = 0
    while next_pts == None and no_searches < max_no_searches:
	next_ptsall, min_est, max_mse_pt = pytave.feval(
            3, "krig_min_find_MADS_oct",
            Am, Jm, curr_bestA, theta, upb, lob,
            N, amin, amax, spc, delta)
	next_pts = check_for_new_points(next_ptsall, Aall)
	no_searches += 1
    return next_pts

def refine(delta, deltamin, spc):
    """Refine the mesh."""
    if delta > deltamin:
	spc = (spc-1.0)*4.0+1.0
        delta = delta/4.0
        return True, delta, spc
    else:
        print "Can not refine anymore with this deltamin"
        return False, None, None

def poll(Aall, Jall, curr_bestA, N, delta, spc, amin, amax):
    """Perform the poll step."""

    poll_pts = None
    try:
	poll_ptsall = pytave.feval(
            1, "MADS_poll_ptsNd3_oct",
            curr_bestA, N, delta, spc, amin, amax)
	poll_ptsall = poll_ptsall[0]
	poll_pts = check_for_new_points(poll_ptsall, Aall)
	poll_pts = pytave.feval(
            1, "sort_poll_pts", Aall, Jall, poll_pts, theta,
            upb, lob, N, amin, amax)

    except Exception as e:
	print " This poll DID NOT work, must refine."

    return poll_pts


def check(Ai, J, nit, curr_bestJ, curr_bestA):
    """Check whether previous step was an improvement."""

    cost_improve = 0
    improve = 0
    Jsurr = []
    Asurr = []
    neval = numpy.zeros((len(J), 1))
    for i in range(0, len(J)):
        if not J[i] >= 1000000:
            Jsurr.append(J[i])
            if Ai.shape ==(3,):
                Asurr.append(Ai)
            else:
                Asurr.append(Ai[i,:])

	if J[i] < curr_bestJ:
            curr_bestJ = J[i]

            if Ai.shape ==(3,):
                curr_bestA = Ai
            else:
                curr_bestA = Ai[i,:]
            curr = [curr_bestJ, curr_bestA, nit, neval, prev_it]
            cost_improve = 1
            improve = 1
	if improve == 0 and nit > 1:
	    cost_improve = 0
	elif nit == 1:
	    cost_improve = 1

    return cost_improve, curr_bestA, curr_bestJ

# add SMF to the octave load path
pytave.feval(0, "addpath", "SMF")

# number of parameters
N = 3

# min and max values for the parameters
amin = numpy.array([100, 0.1, 0.1])
amax = numpy.array([10000, 1.0, 0.5 ])

# initial guess
curr_best = []
curr_bestJ = 1000000
curr_bestA = numpy.array([1000, 0.5, 0.5])

# resolution of the Latin hyper cube
spc = 5*numpy.ones([1,N])
delta = 1.0
deltamin = 1.0/128
ninit = 4*N

# theta and lower and upper bounds on theta
theta = numpy.ones([1,N])*10.0
lob = numpy.ones([1,N])*1e-2
upb = numpy.ones([1,N])*20

neval=ninit
prev_it = 1
cost_improve=1
stop_opt = 0

# number of iterations and max iterations
nit=1
max_nit = 24

Ai = initialize()
Aall, Jall = run_simulations(Ai)

prev_it = "search"
improve = True

converged = False
refine_ok = True

# main loop
while nit <= max_nit and refine_ok and not converged:
   # search step
   if cost_improve:
       Ai_new = search(Aall, Jall, curr_bestA, theta,
                       upb, lob, N, amin, amax, spc, delta)
       prev_it = "search"
       Ai_new = coarsen(Ai_new)
   else:
       # poll step
       if prev_it == "search":
	   Ai_new = poll(Aall, Jall, curr_bestA, N,
                         delta, spc, amin, amax)
	   prev_it = "poll"
       # refine if previous poll did not lead to cost improvement
       if prev_it == "poll":
	   refine_ok, delta, spc = refine(delta, deltamin, spc)
           if refine_ok:
	       Ai_new = search(Aall, Jall, curr_bestA, theta,
                               upb, lob, N, amin, amax, spc, delta)
	       prev_it = "search"
           else:
               Ai_new = None
   nit += 1

   # run simulations on the new parameters
   if not Ai_new == None:
       Ai_new, J_new = run_simulations(Ai_new)

       # stack the new runs to the previous
       Jall = numpy.hstack((Jall, J_new))
       Aall = numpy.vstack((Aall, Ai_new))

       # monitor convergence (write to files)
       monitor(Aall, Jall, nit, curr_bestA, curr_bestJ,
               delta, prev_it, improve, spc)

       # check convergence
       cost_improve, curr_bestA, curr_bestJ = check(
           Ai_new, J_new, nit, curr_bestJ, curr_bestA)
   else:
       cost_improve = 0


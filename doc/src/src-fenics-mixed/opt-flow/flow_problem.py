from dolfin import *

class EssentialBC(Expression):
    def eval(self, v, x):
        if x[0] < DOLFIN_EPS:
            y = x[1]
            v[0] = y*(1-y);  v[1] = 0
        else:
            v[0] = 0; v[1] = 0

    def value_shape(self):
        return (2,)

class K(Expression):
    def __init__(self, K0, x0, c):
        self.K0, self.x0, self.c = K0, x0, c

    def eval(self, v, x):
        x0, K0, c = self.x0, self.K0, self.c
        if abs(x[0] - x0) <= c and abs(x[1] - 0.5) <= c:
            v[0] = K0
        else:
            v[0] = 1

def dirichlet_boundary(x):
    return bool(x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS or \
                x[1] > 1.0 - DOLFIN_EPS)

class FlowProblem2Optimize:
    def __init__(self, K0, x0, c, plot):
        self.K0, self.x0, self.c, self.plot = K0, x0, c, plot

    def run(self):
        K0, x0, c = self.K0, self.x0, self.c

        mesh = UnitSquareMesh(20, 20)
        V = VectorFunctionSpace(mesh, "Lagrange", 2)
        Q = FunctionSpace(mesh, "Lagrange", 1)
        W = MixedFunctionSpace([V,Q])
        u, p = TrialFunctions(W)
        v, q = TestFunctions(W)
        k = K(K0, x0, c)

        u_inflow = EssentialBC()
        bc = DirichletBC(W.sub(0), u_inflow, dirichlet_boundary)
        f = Constant(0)

        a = inner(grad(u), grad(v))*dx + k*inner(u, v)*dx + \
            div(u)*q*dx + div(v)*p*dx
        L = f*q*dx

        w = Function(W)
        solve(a == L, w, bc)

        u, p = w.split()
        u1, u2 = split(u)

        goal1 = assemble(inner(grad(p), grad(p))*dx)
        goal2 = u(1.0, 0.5)[0]*1000
        goal = goal1 + goal2

        if self.plot:
            plot(u)

        key_variables = dict(K0=K0, x0=x0, c=c, goal1=goal1,
                             goal2=goal2, goal=goal)
        print key_variables
        return goal1, goal2

if __name__ == "__main__":
    from sys import argv

    # a test with the porous media placed close to the middle of the domain
    K0=1000; x0 = 0.4; c=0.1

    # near optimal values found by previous SMF run
    # K0 = 5.640625e+02; x0 = 9.156250e-01; c = 1.039063e-01

    p = FlowProblem2Optimize(K0, x0, c, True)
    goal1, goal2 = p.run()

    interactive()


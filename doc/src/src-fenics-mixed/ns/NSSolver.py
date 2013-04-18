from dolfin import *
import math

class NSSolver:
    def __init__(self):
        "initialization"
        mesh = Mesh("aneurysm.xml.gz")
        V = VectorFunctionSpace(mesh, "Lagrange", 2)
        Q = FunctionSpace(mesh, "Lagrange", 1)

        u, p = TrialFunction(V), TrialFunction(Q)
        v, q = TestFunction(V), TestFunction(Q)

        u_1 = Function(V)
        u_ = Function(V)
        p_1 = Function(Q)
        p_ = Function(Q)

        dt = 0.01
        T = 1.0
        nu = 0.0035

        k = Constant(dt)
        f = Constant((0, 0, 0))

        g_noslip = Constant((0, 0, 0))
        bc_noslip = DirichletBC(V, g_noslip, 0)
        self.bcu = [bc_noslip]

        # Form for the velocity step
        F1 = (1/k)*inner(u - u_1, v)*dx + inner(grad(u_1)*u_1, v)*dx + \
             nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
        a1 = lhs(F1)
        self.L1 = rhs(F1)

        # Form for the pressure update
        a2 = inner(grad(p), grad(q))*dx
        self.L2 = -(1/k)*div(u_)*q*dx

        # Form for the velocity update
        a3 = inner(u, v)*dx
        self.L3 = inner(u_, v)*dx - k*inner(grad(p_), v)*dx

        # Assemble matrices
        self.A1 = assemble(a1)
        self.A2 = assemble(a2)
        self.A3 = assemble(a3)

        self.mesh = mesh
        self.u_ = u_
        self.p_ = p_
        self.V = V
        self.Q = Q

        self.ufile = File("results/u.pvd")
        self.pfile = File("results/p.pvd")

    def setIC(self):
        pass

    def inlet_pressure_bc(self, t):
        return math.sin(2*3.14*t)*16000/3 + 16000

    def flux(self):
        n = FacetNormal(self.mesh)
        flux1 = inner(self.u_, n)*ds(1)
        flux2 = inner(self.u_, n)*ds(2)

        Q1 = assemble(flux1)
        Q2 = assemble(flux2)
        Q = [Q1, Q2]
        return Q

    def advance_one_time_step(self, P, t):
        P1, P2 = P
        P0 = 120
        self.t = t

        u_ = self.u_
        p_ = self.p_

        # Compute tentative velocity step
        b1 = assemble(self.L1)
        [bc.apply(self.A1, b1) for bc in self.bcu]
        solve(self.A1, u_.vector(), b1, "gmres", "default")

        # Define pressure boundary conditions
        bc_p_1 = DirichletBC(self.Q, Constant(P[0]), 1)
        bc_p_2 = DirichletBC(self.Q, Constant(P[1]), 2)
        bc_p_3 = DirichletBC(self.Q, self.inlet_pressure_bc(t), 3)
        self.bcp = [bc_p_1, bc_p_2, bc_p_3]

        # Pressure correction
        b2 = assemble(self.L2)
        [bc.apply(self.A2, b2) for bc in self.bcp]
        solve(self.A2, p_.vector(), b2, "gmres", "amg")

        # Velocity correction
        b3 = assemble(self.L3)
        [bc.apply(self.A3, b3) for bc in self.bcu]
        solve(self.A3, u_.vector(), b3, "gmres", "default")

        self.ufile << u_
        self.pfile << p_

        print " max u_ ", u_.vector().max(),
        print " max p_ ", p_.vector().max(),
        print " min p_ ", p_.vector().min()

if __name__ == "__main__":
   solver = NSSolver()
   solver.setIC()
   t = 0
   dt = 0.01
   T = 1.0
   P1, P2 = 0, 0
   while t < T:
       t += dt
       solver.advance_one_time_step((P1, P2), t)


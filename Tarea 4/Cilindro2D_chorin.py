"""
Solucion del problema de Navier-Stokes 2D no estacionario

	-\nu\Delta u + (\nabla u)u+\nabla p = f en ]0,T[\times\Omega
		\nabla\cdot u = 0					en ]0,T[\times\Omega
		   u(0,\cdot) = u_0					en \Omega
					u = u_D					en ]0,T[\times\Gamma_D
		\sigma\cdot n = g					en ]0,T[\times\Gamma_N

utilizando Taylor-Hood P^2-P^1 en el espacio y proyeccion de Chorin


"""

from dolfin import *

# Constantes
mu = 0.01 # Viscocidad del fluido


# Discmretizacion del tiempo
T = 3 # Tiempo final
N = 400 # Pasos
dt = T/N # Tamaño de paso

# Lado derecho
f = Constant((0.0,0.0))

# Cargar malla desde archivo xml(malla hecha con gmsh)
malla = Mesh()
hdf = HDF5File(malla.mpi_comm(), "Mesh/cilindro2D.h5", "r")
hdf.read(malla, "/mesh", False)

dim = 2

dominios = MeshFunction("size_t", malla,dim)
hdf.read(dominios, "/domains")

caras = MeshFunction("size_t", malla, dim-1)
hdf.read(caras, "/bc_markers")

# Definicion espacio de funciones (P^2-P^1)
U = VectorFunctionSpace(malla, "Lagrange", 2) # Velocidad
P = FunctionSpace(malla, "Lagrange", 1) # Presion

# Condiciones de contorno
u_d = Expression(("6*x[1]*(0.41-x[1])/(0.41*0.41)","0.0"), degree=2)
bc_lados = DirichletBC(U,Constant((0,0)), caras, 2)
bc_in = DirichletBC(U, u_d, caras,1)
bc_out = DirichletBC(P, Constant(0), caras, 3)
bcu = [bc_lados, bc_in]
bcp = [bc_out]

n = FacetNormal(malla)
ds = Measure('ds', domain=malla, subdomain_data=caras)

# Funciones test e incognitas
u = TrialFunction(U)
p = TrialFunction(P)
v = TestFunction(U)
q = TestFunction(P)

# Set parameter values
dt = 0.01
T = 3
nu = 0.01

# Define time-dependent pressure boundary condition
#p_in = Expression("sin(3.0*t)", t=0.0)

# Creacion de funciones
u0 = Function(U)
u1 = Function(U)
p0 = Function(P)
p1 = Function(P)

# Define coefficients
k = Constant(dt)

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Condiciones de contorno
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Creacion VTKs
ufile = File("Salidas_Cilindro2D_chorin/velocidad.pvd")
pfile = File("Salidas_Cilindro2D_chorin/presion.pvd")

u1.rename('Velocidad','Velocidad del fluido')
p1.rename('Presion','Presion del fluido')

# Precondicionador
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    #p_in.t = t

    # Aproximacion de la velocidad
	b1 = assemble(L1)
	[bc.apply(b1) for bc in bcu]
	solve(A1, u1.vector(), b1, 'bicgstab', 'hypre_amg')

	# Calculo de la presion
	b2 = assemble(L2)
	[bc.apply(b2) for bc in bcp]
	solve(A2, p1.vector(), b2, 'bicgstab', 'hypre_amg')

	# Calculo de la velocidad
	b3 = assemble(L3) 
	solve(A3, u1.vector(), b3, 'cg' , 'sor' )
	
	ufile << u1
	pfile << p1

	u0.assign(u1)
	p0.assign(p1)
	t += dt

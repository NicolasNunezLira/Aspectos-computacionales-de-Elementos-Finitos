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
from numpy import pi

# Constantes
mu = 0.01 # Viscocidad del fluido

# Discmretizacion del tiempo
T = 100 # Tiempo final
N = 3000 # Pasos
dt = T/N # Tamaño de paso

# Lado derecho
f = Constant((0.0,0.0))
#g = Constant((1,0))

# Cargar malla desde archivo xml(malla hecha con gmsh)
malla = Mesh()
hdf = HDF5File(malla.mpi_comm(), "Mesh/Square.h5", "r")
hdf.read(malla, "/mesh", False)

dim = 2

dominios = MeshFunction("size_t", malla,dim)
hdf.read(dominios, "/domains")

caras = MeshFunction("size_t", malla, dim-1)
hdf.read(caras, "/bc_markers")


# Definicion de espacios (P^2-P^1)
U = VectorFunctionSpace(malla,'Lagrange',2) #Velocidad
P = FunctionSpace(malla,'Lagrange',1) #Presion

# Condiciones de contorno
u_d = Expression(("cos(pi*x[0])*sin(pi*x[1])*exp(-2*t*mu*pow(pi,2))","cos(pi*x[1])*sin(pi*x[0])*exp(-2*t*mu*pow(pi,2))"),mu=mu,pi=pi,t=0, degree=2)
p_d = Expression("-0.25*(cos(2*pi*x[0])+cos(2*pi*x[1]))*exp(-4*t*mu*pow(pi,2))",mu=mu,pi=pi,t=0, degree=1)
bc_lados = DirichletBC(U,u_d, caras ,1)
bc_out = DirichletBC(P, p_d, caras, 1)
bcu = [bc_lados]
bcp = [bc_out]

n = FacetNormal(malla)
ds = Measure('ds', domain=malla, subdomain_data=caras)

#Funciones test e incognitas
u = TrialFunction(U)
p = TrialFunction(P)
v = TestFunction(U)
q = TestFunction(P)

# Creacion de las funciones
u0 = Function(U)
u1 = Function(U)
p0 = Function(P)
p1 = Function(P)

# Define coefficients
k = Constant(dt)

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     mu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
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
ufile = File("Salidas_Square_chorin/velocidad.pvd")
pfile = File("Salidas_Square_chorin/presion.pvd")

u1.rename('Velocidad','Velocidad del fluido')
p1.rename('Presion','Presion del fluido')

# Precondicionador
#prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
t = 0
for i in range(N+1):
	print(t)
	# Aproximacion de la velocidad
	A1 = assemble(a1)
	[bc.apply(A1) for bc in bcu]
	b1 = assemble(L1)
	[bc.apply(b1) for bc in bcu]
	solve(A1, u1.vector(), b1, 'bicgstab', 'hypre_amg')

	# Calculo de la presion
	A1 = assemble(a2)
	[bc.apply(A2) for bc in bcp]
	b2 = assemble(L2)
	[bc.apply(b2) for bc in bcp]
	solve(A2, p1.vector(), b2, 'bicgstab', 'hypre_amg')

	# Calculo de la velocidad
	A3 = assemble(a3)
	b3 = assemble(L3) 
	solve(A3, u1.vector(), b3, 'cg' , 'sor' )
	
	ufile << u1
	pfile << p1

	u0.assign(u1)
	p0.assign(p1)
	t += dt
	u_d.t = t
	p_d.t = t

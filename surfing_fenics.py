###########################################################################################################################################################################################################################
# This FEniCS code implements the phase-field model of Kumar et al., JMPS 2020 to solve the 'Surfing problem' to show that fracture propagation follows the Griffith principle.
#
# Input: Set the material properties in lines 17-23. Create a geometry with a crack substantially bigger than Irwin's characteristic length scale. 
#		 Set appropriate surfing boundary conditions in line 68.
#
# Output: Energy release rate, G, is calculated by computing J-integral in a loop around the crack. It is printed in the file 'Surfing.txt'. G is expected to stabilize around Gc.
#		  The crack growth rate can be calculated from visualizing the results in paraview.
#
# Contact Aditya Kumar (akumar355@gatech.edu) for questions.
###########################################################################################################################################################################################################################


from dolfin import *
import numpy as np
import time
import sys

# Material properties
E, nu = 9800, 0.13
mu, kappa = E/(2*(1 + nu)), E/(3*(1 - 2*nu))
lmbda = kappa - 2 / 3 * mu
Gc= 0.091125	
sts, scs= 27, 77	

lch=3*Gc*E/8/(sts**2)

eps=0.35      

shs = (2/3)*sts*scs/(scs-sts)
Wts = sts**2/(2*E)
Whs = shs**2/(2*kappa)


comm = MPI.comm_world 
comm_rank = MPI.rank(comm)

# Create mesh and define function space
mesh=Mesh("mesh/surfing_mesh.xml") 
h=FacetArea(mesh)          #area/length of a cell facet on a given mesh
h_avg = (h('+') + h('-'))/2
n=FacetNormal(mesh)

delta = (1+3*h/(8*eps))**(-2) * ((sts + (1+2*np.sqrt(3))*shs)/((8+3*np.sqrt(3))*shs)) * 3*Gc/(16*Wts*eps) + (1+3*h/(8*eps))**(-1) * (2/5)

# Choose phase-field model
phase_model=1;  #1 for linear model, 2 for quadratic model

V = VectorFunctionSpace(mesh, "CG", 1)   #Function space for u
Y = FunctionSpace(mesh, "CG", 1)         #Function space for z


# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side, tol) && abs(x[1]) > tol && on_boundary", side = 0.0, tol=1e-4)
right =  CompiledSubDomain("near(x[0], side, tol) && on_boundary", side = 30.0, tol=1e-4)  #30.0
bottom =  CompiledSubDomain("near(x[1], side, tol) && on_boundary", side = -5.0, tol=1e-4)  #5.0
top =  CompiledSubDomain("near(x[1], side, tol) && on_boundary", side = 5.0, tol=1e-4)
righttop = CompiledSubDomain("abs(x[0]-30.0)<1e-4 && abs(x[1]-5.0)<1e-4 ")
corner = CompiledSubDomain("abs(x[1]-0.0)<1e-4 && x[0]<5+1e-4 && x[0]>4.5")

#############################################################	
set_log_level(40)  #Error level=40, warning level=30
parameters["linear_algebra_backend"] = "PETSc"
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}	
############################################################
	
##################################################################################
# Define Dirichlet boundary conditions
##################################################################################
c=Expression("K1/(2*mu)*sqrt(sqrt(pow(x[0]-V*(t+0.1),2)+pow(x[1],2))/(2*pi))*(kap-cos(atan2(x[1],(x[0]-V*(t+0.1)))))*cos(atan2(x[1],(x[0]-V*(t+0.1)))/2)",degree=4,t=0, V=20, K1=30, mu=4336.28, kap=2.54)
r=Expression("K1/(2*mu)*sqrt(sqrt(pow(x[0]-V*(t+0.1),2)+pow(x[1],2))/(2*pi))*(kap-cos(atan2(x[1],(x[0]-V*(t+0.1)))))*sin(atan2(x[1],(x[0]-V*(t+0.1)))/2)",degree=4,t=0, V=20, K1=30, mu=4336.28, kap=2.54)

bcl= DirichletBC(V.sub(0), Constant(0.0), righttop, method='pointwise'  )
bcb2 = DirichletBC(V.sub(1), r, bottom )
bct2 = DirichletBC(V.sub(1), r, top)
bcs = [ bcl, bcb2,  bct2]

cz=Constant(1.0)
bcb_z = DirichletBC(Y, cz, bottom)
bct_z = DirichletBC(Y, cz, top)
cz2=Constant(0.0)
bcc_z = DirichletBC(Y, cz2, corner)
bcs_z=[ bcb_z, bct_z, bcc_z]

########################################################################
# Define functions
########################################################################
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
u_inc = Function(V)
dz = TrialFunction(Y)            # Incremental phase field
y  = TestFunction(Y)             # Test function
z  = Function(Y)                 # Phase field from previous iteration
z_check = Function(Y)
z_inc = Function(Y)
d = u.geometric_dimension()

##############################################################
#Initialisation of displacement field,u and the phase field,z
##############################################################
u_init = Constant((0.0,  0.0))
u.interpolate(u_init)
for bc in bcs:
	bc.apply(u.vector())

z_init = Constant(1.0)
z.interpolate(z_init)
for bc in bcs_z:
	bc.apply(z.vector())

z_ub = Function(Y)
z_ub.interpolate(Constant(1.0))	
z_lb = Function(Y)
z_lb.interpolate(Constant(-0.0))
	
u_prev = Function(V)
assign(u_prev,u)
z_prev = Function(Y)
assign(z_prev,z)
	
#################################################
###Label the dofs on boundary
#################################################
def extract_dofs_boundary(V, bsubd):	
	label = Function(V)
	label_bc_bsubd = DirichletBC(V, Constant((1,1)), bsubd)
	label_bc_bsubd.apply(label.vector())
	bsubd_dofs = np.where(label.vector()==1)[0]
	return bsubd_dofs

#Dofs on which reaction is calculated
top_dofs=extract_dofs_boundary(V,top)
y_dofs_top=top_dofs[1::d]


boundary_subdomains = MeshFunction("size_t", mesh, 1)
boundary_subdomains.set_all(0)
left.mark(boundary_subdomains,1)
right.mark(boundary_subdomains,1)
bottom.mark(boundary_subdomains,2)
top.mark(boundary_subdomains,2)


# Define new measures associated with the interior domains 	
ds = ds(subdomain_data=boundary_subdomains)


# Elasticity parameters
mu, lmbda, kappa = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu))), Constant(E/(3*(1 - 2*nu)))

def local_project(v, V, u=None):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return
	
def energy(v):
	return mu*(inner(sym(grad(v)),sym(grad(v))) + ((nu/(1-nu))**2)*(tr(sym(grad(v))))**2 )+  0.5*(lmbda)*(tr(sym(grad(v)))*(1-2*nu)/(1-nu))**2 
	
def epsilon(v):
	return sym(grad(v))

def sigma(v):
	return 2.0*mu*sym(grad(v)) + (lmbda)*tr(sym(grad(v)))*(1-2*nu)/(1-nu)*Identity(len(v))

def sigmavm(v):
	sig = sigma(v)
	return sqrt(1/2*(inner(sig-1/3*tr(sig)*Identity(len(v)), sig-1/3*tr(sig)*Identity(len(v))) + (1/9)*tr(sig)**2 ))

##configurational external force to model strength
beta1 = -(delta*Gc)/(shs*8*eps) + (2*Whs)/(3*shs)
beta2 = -(3**0.5*(3*shs - sts)*delta*Gc)/(shs*sts*8*eps) - (2*Whs)/(3**0.5*shs) + (2*3**0.5*Wts)/(sts)

pen=1000*(3*Gc/8/eps) * conditional(lt(delta, 1), 1, delta)
eta=1e-5	

# Stored strain energy density (compressible L-P model)
psi1 =(z**2+eta)*(energy(u))	
psi11=energy(u)
# stress=(z**2+eta)*sigma(u)

I1_d = (z**2+eta)*tr(sigma(u))
SQJ2_d = (z**2+eta)*sigmavm(u)
ce = beta2*SQJ2_d + beta1*I1_d + z*(1-sqrt(I1_d**2)/I1_d)*psi11 

# Total potential energy
Pi = psi1*dx 

# Compute first variation of Pi (directional derivative about u in the direction of v)
R = derivative(Pi, u, v)

# Compute Jacobian of R
Jac = derivative(R, u, du) 

#To use later for memory allocation for these tensors
A=PETScMatrix()
b=PETScVector()

#Balance of configurational forces PDE
Wv=pen/2*((abs(z)-z)**2 + (abs(1-z) - (1-z))**2 )*dx
Wv2=conditional(le(z, 0.05), 1, 0)*40*pen/2*( 1/4*( abs(z_prev-z)-(z_prev-z) )**2 )*dx

if phase_model==1:
	R_z = y*2*z*(psi11)*dx - y*(ce)*dx +3*delta*Gc/8*(y*(-1)/eps + 2*eps*inner(grad(z),grad(y)))*dx + derivative(Wv,z,y) #+  derivative(Wv2,z,y)  #linear model
else:
	R_z = y*2*z*(psi11)*dx+ y*(ce)*dx+ Gc*(y*(z-1)/eps + eps*inner(grad(z),grad(y)))*dx + derivative(Wv2,z,y)  #quadratic model
	
# Compute Jacobian of R_z
Jac_z = derivative(R_z, z, dz)



## Define the solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "cg",   
                                          "preconditioner": "amg",						  
                                          "maximum_iterations": 10,
                                          "report": True,
                                          "error_on_nonconvergence": False}}			


#time-stepping parameters
T=1

Totalsteps=50
startstepsize=1/Totalsteps
stepsize=startstepsize
t=stepsize
step=1
rtol=1e-8
printsteps=5 
tau=0

start_time=time.time()
while t-stepsize < T:

	if comm_rank==0:
		print('\n \nStep= %d' %step, 't= %f' %t, 'Stepsize= %e' %stepsize)
		sys.stdout.flush()
	
	c.t=t; r.t=t
	stag_iter=1
	rnorm_stag=1
	while stag_iter<50 and rnorm_stag > rtol:
		start_time=time.time()
		##############################################################
		#First PDE
		##############################################################		
		Problem_u = NonlinearVariationalProblem(R, u, bcs, J=Jac)
		solver_u  = NonlinearVariationalSolver(Problem_u)
		solver_u.parameters.update(snes_solver_parameters)
		(iter, converged) = solver_u.solve()
		
		##############################################################
		#Second PDE
		##############################################################
		assign(z_check,z)
		Problem_z = NonlinearVariationalProblem(R_z, z, bcs_z, J=Jac_z)
		solver_z  = NonlinearVariationalSolver(Problem_z)
		solver_z.parameters.update(snes_solver_parameters)
		(iter, converged) = solver_z.solve()
		##############################################################
		
		
		if comm_rank==0:
			print("--- %s seconds ---" % (time.time() - start_time))
			sys.stdout.flush()

	  
		###############################################################
		#Residual check for stag loop
		###############################################################
		# b=assemble(-R, tensor=b)
		# for bc in bcs:
		# 	bc.apply(b)
		# rnorm_stag=b.norm('l2')	
		rnorm_stag = norm(z_check.vector() - z.vector())
		if comm_rank==0:
			print('Stag Iteration no= %d' %stag_iter,  'Residual= %e' %rnorm_stag)
			sys.stdout.flush()
		stag_iter+=1  

	
	######################################################################
	#Post-Processing
	######################################################################
	assign(u_prev,u)
	assign(z_prev,z)
	
	tau+=stepsize

		
	####Calculate Reaction
	stress=(z**2+eta)*sigma(u)
	JI1=(psi1-dot(dot(stress,n),u.dx(0)))*ds(1)
	JI2=(-dot(dot(stress,n),u.dx(0)))*ds(2)
	Jintegral=assemble(JI1)+assemble(JI2)

	if comm_rank==0:
		with open('Surfing_Jint.txt', 'a') as rfile:
			rfile.write("%s %s\n" % (str(t), str(Jintegral)))
		

	####Plot solution on incremental steps
	if step % 1==0:
		file_results = XDMFFile( "Paraview/Surfing_step_" + str(step) + ".xdmf" )
		file_results.parameters["flush_output"] = True
		file_results.parameters["functions_share_mesh"] = True
		u.rename("u", "displacement field")
		z.rename("z", "phase field")
		file_results.write(u,t)
		file_results.write(z,t)
	
	#time stepping
	step+=1
	t+=stepsize
 
#######################################################end of all loops
if comm_rank==0:	
	print("--- %s seconds ---" % (time.time() - start_time))


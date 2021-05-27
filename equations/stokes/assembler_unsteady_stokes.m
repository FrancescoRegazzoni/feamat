function [H,B,b] = assembler_unsteady_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions,neumann_functions,dt)
% Assemble steady unsteady stokes matrix and rhs with boundary conditions.
% Backward Euler method is used for time discretization.
% input=
%           fespace_u: finite elemnet space for velocity
%           fespace_p: finite elemnet space for pressure
%           fun: anonymous function or scalar of the forcing term
%               If scalar, fun must be [0;0] and the code is optimized
%           nu: anonymous function or scalar for the diffusion term. If
%           scalar, the code is optimized on structured meshes
%           dirichlet_functions: Dirichlet boundary data
%           neumann_functions: Neumann_boundary data
%           dt: Time step size
% output=
%           H: system matrix
%           B: right-hand side matrix (multiplying state at previous time)
%           b: right-hand side

[H,b] = assembler_steady_stokes(fespace_u,fespace_p,fun,nu,dirichlet_functions,neumann_functions);

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);
n_tot = 2 * n_nodes_u + n_nodes_p;
indices_u1 = 1:n_nodes_u;
indices_u2 = n_nodes_u+1:n_nodes_u*2;
B = sparse(n_tot,n_tot);
M = assemble_mass(fespace_u);
M = apply_dirichlet_bc_matrix(M,fespace_u,0);
B(indices_u1,indices_u1) = M/dt;
B(indices_u2,indices_u2) = M/dt;

H = B + H;
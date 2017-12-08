function [A,b] = assembler_poisson(fespace,fun,mu,dirichlet_functions,neumann_functions)
% Assemble poisson matrix with boundary conditions
% input=
%           fespace: finite elemnet space
%           fun: anonymous function of the forcing term
%           mu: anonymous function or scalar of the diffusion coefficient
%               If scalar, the code is optimized on structured meshes
%           dirichlet_functions: Dirichlet boundary data
%           neumann_functions: Neumann_boundary data
% output=
%           A: system matrix
%           b: right handside

bc_flags = fespace.bc;

thereisneumann = 1;

if (length(find(bc_flags)) == 4)
    thereisneumann = 0;
end

A = assemble_stiffness(mu,fespace);
b = assemble_rhs(fespace,fun);

if (thereisneumann)
   b = apply_neumann_bc(b,fespace,neumann_functions); 
end

% Apply Dirichlet boundary conditions
[A,b] = apply_dirichlet_bc(A,b,fespace,dirichlet_functions);



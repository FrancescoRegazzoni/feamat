clear all
clc
fprintf('initializeing example...\n')

bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1;
H = 1;

n_elements_x = 20;
n_elements_y = 20;

mesh = create_mesh(bottom_left_corner_x, ...
                   bottom_left_corner_y, ...
                   L,H,n_elements_x,n_elements_y);

bc_flags = [1 1 1 1];

fespace_u = create_fespace(mesh,'P2',bc_flags);
fespace_p = create_fespace(mesh,'P1',bc_flags);

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);
n_tot = 2 * n_nodes_u + n_nodes_p;

f = [0;0];
mu = 1;

U = 500;

dt = 1e-3;
T = 0.1;

dirichlet_functions = @(x) [0 0;0 0;U 0;0 0]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

fprintf('assembling matrices... ')
[A,B,b] = assembler_unsteady_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions,neumann_functions,dt);
fprintf('done!\n')

tt = 0:dt:T;

figure();

vecsol = zeros(n_tot,1);

for iT = 2:length(tt)
    fprintf('iter %d\n',iT)
    fprintf('   Solving linear system... ')
    
    u = vecsol(1:2*n_nodes_u);
    A_C = add_convective_term(A,u,fespace_u);
    vecsol = A_C \ (B*vecsol+ b);
    
    fprintf('done!\n')

    sol = get_fluid_solution(vecsol,fespace_u,fespace_p);
    plot_fe_fluid_function(sol,'U');
    pause(1e-16);
    export_vtk_fluid(sol,sprintf('example_unsteady_navier_stokes_%06d.vtk', iT))
end
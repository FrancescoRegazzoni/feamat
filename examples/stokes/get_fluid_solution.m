function sol = get_fluid_sol(vecsol,fespace_u,fespace_p)

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);

indices_u1 = 1:n_nodes_u;
indices_u2 = n_nodes_u+1:2*n_nodes_u;
indices_p = 2*n_nodes_u+1:2*n_nodes_u+n_nodes_p;

sol.n_nodes_u = n_nodes_u;
sol.n_nodes_p = n_nodes_p;
sol.u1 = vecsol(indices_u1);
sol.u2 = vecsol(indices_u2);
sol.p = vecsol(indices_p);
sol.fespace_u = fespace_u;
sol.fespace_p = fespace_p;
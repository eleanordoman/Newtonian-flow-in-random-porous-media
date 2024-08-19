function [newtonian_solve, network_periodic] = newtonian_flow_solver_impose_pressure_drop_periodic(network, newtonian_solve, domain, dP, mu)
%Solve for pressures at nodes/pores and fluxes along edges using resistance input into function. A is incidence
%matrix of system. System is periodic.

%define peridoic incidence matrix ie

k_top=0; top_vertices=[];
k_bottom = 0; bottom_vertices=[];
for i=1:network.K
    if network.pores(i, 2)==0
        k_bottom=k_bottom+1;
        bottom_vertices = [bottom_vertices; i];
    elseif network.pores(i, 2)==domain.size_y
        k_top=k_top+1;
        top_vertices=[top_vertices; i];
    end
end

ind1 = find(network.edgeY(:, 1)==0);
ind2 = find(network.edgeY(:, 2)==0);
bottom_edges = intersect(ind1 ,ind2); m_bottom = length(bottom_edges);
ind3 = find(network.edgeY(:, 1)==domain.size_y);
ind4 = find(network.edgeY(:, 2)==domain.size_y);
top_edges = intersect(ind3, ind4); m_top = length(top_edges);

network_periodic.pores = network.pores;
network_periodic.pores(top_vertices, :)=[];

network_periodic.edge_width = network.edge_width;
edgeX_top= network.edgeX(top_edges, :);
top_vertices_to_remove=[];
for i= bottom_edges.'
    indx = find(network.edgeX(i, 1)== edgeX_top(:, 1));
    indy = find(network.edgeX(i, 2)== edgeX_top(:, 2));
    ind = intersect(indx, indy);
    if isempty(ind)==1
        indx = find(network.edgeX(i, 1)== edgeX_top(:, 2));
        indy = find(network.edgeX(i, 2)== edgeX_top(:, 1));
        ind = intersect(indx, indy);
    end
    if isempty(ind)==0
        network_periodic.edge_width(i) = network.edge_width(i)+ network.edge_width(top_edges(ind));
        top_vertices_to_remove=[top_vertices_to_remove, top_edges(ind)];
    end
end

for i= bottom_edges.'
    newtonian_solve.R(i) = (mu*pi^4*domain.r^0.5/(8*domain.depth^2*network_periodic.edge_width(i)^1.5))*g(network_periodic.edge_width(i)/domain.depth, network.edge_length(i)/(network_periodic.edge_width(i)*domain.r)^0.5);
end
newtonian_solve.R(top_vertices_to_remove)=[];
network_periodic.edge_width(top_vertices_to_remove) = [];

network_periodic.edgeX = network.edgeX;
network_periodic.edgeY = network.edgeY;
network_periodic.edgeX(top_vertices_to_remove, :) = [];
network_periodic.edgeY(top_vertices_to_remove, :) = [];

network_periodic.edge_length = network.edge_length;
network_periodic.edge_length(top_vertices_to_remove) = [];

network_periodic.M = length(network_periodic.edge_width);
network_periodic.m = network.m; network_periodic.m1 = network.m1; network_periodic.m2 = network.m2;

network_periodic.K = length(network_periodic.pores);
network_periodic.k = network.k; network_periodic.k1 = network.k1; network_periodic.k2 = network.k2;

A=zeros(size(network_periodic.edgeX, 1), length(network_periodic.pores));
for i=1:size(network_periodic.edgeX, 1) %determine incidence matrix
    for j=1:length(network_periodic.pores)
        if round(network_periodic.edgeX(i, 1), 4)==round(network_periodic.pores(j, 1), 4)
            if round(network_periodic.edgeY(i, 1), 4)==round(network_periodic.pores(j, 2), 4)
                A(i, j)=-1;
            elseif round(network_periodic.edgeX(i, 1), 4) == round(network_periodic.edgeX(i, 2), 4)
                if round(network_periodic.edgeY(i, 2), 4)==round(network_periodic.pores(j, 2), 4)
                    A(i, j)=1;
                elseif network_periodic.edgeY(i, 1) == 0 && round(network_periodic.pores(j, 2), 4)==domain.size_y
                    A(i, j)=-1;
                elseif network_periodic.edgeY(i, 2) == 0 && round(network_periodic.pores(j, 2), 4)==domain.size_y
                    A(i, j)=1;
                elseif round(network_periodic.edgeY(i, 1), 4) == domain.size_y && network_periodic.pores(j, 2)== 0
                    A(i, j)=-1;
                elseif round(network_periodic.edgeY(i, 2), 4) == domain.size_y && network_periodic.pores(j, 2)== 0
                    A(i, j)=1;
                end
            elseif network_periodic.edgeY(i, 1) == 0 && round(network_periodic.pores(j, 2), 4)==domain.size_y
                A(i, j)=-1;
            elseif round(network_periodic.edgeY(i, 1), 4) == domain.size_y && network_periodic.pores(j, 2)== 0
                A(i, j)=-1;
            end
        elseif round(network_periodic.edgeX(i, 2), 4)==round(network_periodic.pores(j, 1), 4)
            if round(network_periodic.edgeY(i, 2), 4)==round(network_periodic.pores(j, 2), 4)
                A(i, j)=1;
            elseif network_periodic.edgeY(i, 2) ==0 && round(network_periodic.pores(j, 2), 4)== domain.size_y
                A(i, j)=1;
            elseif round(network_periodic.edgeY(i, 2), 4) ==domain.size_y && network_periodic.pores(j, 2)== 0
                A(i, j)=1;
            end
        end
    end
end

network_periodic.A = A;

A_e = network_periodic.A(1:network_periodic.m, 1:network_periodic.k);
A_i = network_periodic.A(network_periodic.m+1:end, network_periodic.k+1:end);
Abar = network_periodic.A(1:network_periodic.m, network_periodic.k+1:end);

R_e=diag(newtonian_solve.R(1:network_periodic.m));
R_i=diag(newtonian_solve.R(network_periodic.m+1:end));

P_e=dP*[ones(network_periodic.k1, 1); zeros(network_periodic.k2, 1)];

C = [Abar, R_e, zeros(network_periodic.m, network_periodic.M-network_periodic.m);
    A_i, zeros(network_periodic.M-network_periodic.m, network_periodic.m), R_i;
    zeros(network_periodic.K-network_periodic.k, network_periodic.K-network_periodic.k), Abar.', A_i.'];

newtonian_solve.C = C;

X = C\[-A_e*P_e; zeros(network_periodic.M+network_periodic.K-network_periodic.k-network_periodic.m, 1)];

newtonian_solve.P = [P_e; X(1:network_periodic.K-network_periodic.k)];
newtonian_solve.Q = [X(network_periodic.K-network_periodic.k+1:end)];
newtonian_solve.dP = dP;

end

function f = f(x)
% x is ratio between width and height of pipe ie x=h/D
    n=0:500;
    f_holder=[];
    for i=1:max(n)+1
        f_holder = [f_holder; (x-2*tanh((2*n(i)+1)*x*0.5*pi)./((2*n(i)+1)*pi))./((2*n(i)+1)^4*x.^2)];
    end
    f = sum(f_holder);
end

function g = g(x, y)
%x is ratio between minium gap and height of pipe ie x = h/D
    fun = @(xi) 1./((1+xi.^2).^2.*f(x.*(1+xi.^2)));
    g = integral(fun, -y/2, y/2);
end


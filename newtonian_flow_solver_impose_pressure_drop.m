function newtonian_solve = newtonian_flow_solver_impose_pressure_drop(network, newtonian_solve, dP)
%Solve for pressures at nodes/pores and fluxes along edges using resistance input into function. A is incidence
%matrix of system. System is not periodic, with no-flow imposed imposed at
%top and botom of domain.

K=network.K; k=network.k; M=network.M; m=network.m;
k1=network.k1; k2=network.k2;

A_e = network.A(1:m, 1:k);
A_i = network.A(m+1:end, k+1:end);
Abar = network.A(1:m, k+1:end);

R_e=diag(newtonian_solve.R(1:m));
R_i=diag(newtonian_solve.R(m+1:M));

P_e=dP*[ones(k1, 1); zeros(k2, 1)];

C = [Abar, R_e, zeros(m, M-m);
    A_i, zeros(M-m, m), R_i;
    zeros(K-k, K-k), Abar.', A_i.'];

newtonian_solve.C = C;

X = C\[-A_e*P_e; zeros(M+K-k-m, 1)];

newtonian_solve.P = [P_e; X(1:K-k)];
newtonian_solve.Q = [X(K-k+1:end)];
newtonian_solve.dP = dP;

Qin = sum(newtonian_solve.Q(1:network.k1));
Qout = sum(newtonian_solve.Q(network.k1+1:network.k));

if round(abs(Qin-Qout), 4)~=0
    disp('Error: mass is not conserved')
end

end


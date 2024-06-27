%generate disordered geometries 
%Note: units are arbitary until Newtonian simulation is calculated
clear all;

domain.r = 5;  %radius of inclusions (um)
domain.size_x = 300; %length of array (um)
domain.size_y=100; %width of array (um)
domain.min_width=2; %minium  distances between two inclusions (um)
domain.depth = 6; %domain depth (um)

%to allow inclusions to overlap with boundary of domain set
domain.min_boundary = -1;

%otherwise to impose minium distance between inclusion and boundary of
%array uncomment following:
% domain.min_boundary = 5; %minium distance between inclusion and boundary
% of array (um); if min_boundary<0 inclusions can overlap with edges

domain.pixel_scaling = 2; %number of pixels per um - increase for better resolution however will also increase compuational time

domain = monodisperse_generator(domain); %generate geometry

%plot domain
figure()
spy(domain.image)

%evaluate the porosity based on pixel image
porosity = 1-sum(sum(domain.image))/(size(domain.image, 1)*size(domain.image, 2)); 

%to map the network:
% to allow vertical lines: vertical_lines = 0, else vertical_lines = 1
vertical_lines=0;
network = map_network(domain, vertical_lines);

%plot pores and edges:
figure()
spy(domain.image)
hold on
plot(network.pores(:, 1)*domain.pixel_scaling, network.pores(:, 2)*domain.pixel_scaling, 'o')
for i=1:network.M
    plot(network.edgeX(i, :)*domain.pixel_scaling, network.edgeY(i, :)*domain.pixel_scaling)
end

%Newtonian simulation:

%viscosity & pressure drop over domain:
mu = 1; %Pa s
dP = 500; %Pa

%calculate vessel resistance:
newtonian_solve = calculate_vessel_resistance(domain, network, mu);

%calculate internal pressures & fluxes throughout domain for SOLID
%boundaries at top and bottom of domain
newtonian_solve = newtonian_flow_solver_impose_pressure_drop(network, newtonian_solve, dP); 

%check total flux in is same as total flux out:
Qin = sum(abs(newtonian_solve.Q(1:network.k1)));
Qout = sum(abs(newtonian_solve.Q(network.k1+1:network.k)));

%plot pressures and fluxes
plot_flow_and_pressure(domain, network, newtonian_solve)
%ensure you expand figure windows or colour bar can display wrong max value

%calculate internal pressures & fluxes throughout domain for PERIODIC
%boundaries at top and bottom of domain
[newtonian_solve_periodic, network_periodic] = newtonian_flow_solver_impose_pressure_drop_periodic(network, newtonian_solve, domain, dP, mu);

%check total flux in is same as total flux out:
Qin_periodic = sum(abs(newtonian_solve_periodic.Q(1:network_periodic.k1)));
Qout_periodic = sum(abs(newtonian_solve_periodic.Q(network.k1+1:network_periodic.k)));

%plot pressures and fluxes for periodic case
plot_flow_and_pressure(domain, network_periodic, newtonian_solve_periodic)


function plot_flow_and_pressure(domain, network, newtonian_solve)
%Plot flow paths and pressures for flow_solver_lubrication_theory

edgeX = network.edgeX*domain.pixel_scaling;
edgeY = network.edgeY*domain.pixel_scaling;
pores = network.pores*domain.pixel_scaling;

Q = newtonian_solve.Q;
P = newtonian_solve.P;

scaling_factor = domain.pixel_scaling;

%plot pressure plot
X = domain.X*scaling_factor;
Y = domain.Y*scaling_factor;

imagesizex=domain.size_x*scaling_factor;
imagesizey=domain.size_y*scaling_factor;

expanded_X=[ X-imagesizex, X-imagesizex, X-imagesizex, X, X, X, X+imagesizex, X+imagesizex, X+imagesizex];
expanded_Y=[Y-imagesizey, Y, Y+imagesizey, Y-imagesizey, Y, Y+imagesizey, Y-imagesizey, Y, Y+imagesizey];
positions1=unique([expanded_X.', expanded_Y.'], 'rows'); %remove dublicate cylinder positions
expanded_X=positions1(:, 1).';
expanded_Y=positions1(:, 2).';

dt = delaunayTriangulation(expanded_X.', expanded_Y.');
face_centres = circumcenter(dt);
Alpha=zeros(size(face_centres, 1), network.K);
problem_vertices=[];
faces_picked = [];
for i = 1:network.K
    indx = find(round(pores(i, 1)/100)==round(face_centres(:, 1)/100));
    indy = find(round(pores(i, 2), 3)==round(face_centres(:, 2), 3));
    ind = intersect(indx, indy);
    if length(ind)==1
        Alpha(ind, i)=1;
        faces_picked=[faces_picked, ind];
    elseif length(ind)>1
        trial_faces = face_centres(ind, :);
        [xa, ya] = closestpoint(pores(i, :), trial_faces);
        indy = find(face_centres(:, 2)==ya);
        if length(indy)>1
            indx = find(face_centres(:, 1)==xa);
            ind = intersect(indx, indy);
        else
            ind = indy;
        end
        Alpha(ind, i)=1;
        faces_picked=[faces_picked, ind];
    elseif isempty(ind)
        problem_vertices = [problem_vertices, i];
    end
end

unclaimed_faces = face_centres;
unclaimed_faces(faces_picked, :) = [];
for i=problem_vertices
    [xf_point, yf_point ] = closestpoint(pores(i, :), unclaimed_faces);
    indy = find(face_centres(:, 2)==yf_point);
    if length(indy)>1
        indx = find(face_centres(:, 1)==xf_point);
        ind = intersect(indx, indy);
    else
        ind = indy;
    end
    Alpha(ind, i)=1;
end

max_P = max(P);
min_P = min(P);
range_P = max_P-min_P;

figure()
axis off
hold on;

for i=1:network.K
    ind = find(Alpha(:, i)==1);
    if isempty(ind)==0
        points = dt.ConnectivityList(ind, :);
        x1 = dt.Points(points(1), 1); y1 = dt.Points(points(1), 2);
        x2 = dt.Points(points(2), 1); y2 = dt.Points(points(2), 2);
        x3 = dt.Points(points(3), 1); y3 = dt.Points(points(3), 2);
        
        patch([x1; x2; x3], [y1; y2; y3], [0, (max_P-P(i))/range_P, 1-(max_P-P(i))/(2*range_P)])
    end
end
colormap(flipud(winter))
a = colorbar('XTickLabel', {num2str(min_P), num2str(min_P +0.1*range_P), num2str(min_P +0.2*range_P), num2str(min_P +0.3*range_P), num2str(min_P +0.4*range_P), num2str(min_P +0.5*range_P), num2str(min_P +0.6*range_P), num2str(min_P +0.7*range_P), num2str(min_P +0.8*range_P), num2str(min_P +0.9*range_P), num2str(min_P +range_P)});
ylabel(a, 'Pressure, Pa')

spy(domain.image, 'k')

%plot fluxes

figure()
spy(domain.image, 'k')
axis off
hold on;

max_Q=max(abs(Q));
for i=1:network.M
    plot(edgeX(i, :), edgeY(i, :), 'color', [1, 1-abs(Q(i))/max_Q, 0], 'linewidth', 5)
end
colormap(flipud(autumn))
a = colorbar('XTickLabel', {num2str(0), num2str(0.1*max_Q), num2str(0.2*max_Q), num2str(0.3*max_Q), num2str(0.4*max_Q), num2str(0.5*max_Q), num2str(0.6*max_Q), num2str(0.7*max_Q), num2str(0.8*max_Q), num2str(0.9*max_Q), num2str(max_Q)});
ylabel(a, 'Volume Flux, um^3/s ')

end
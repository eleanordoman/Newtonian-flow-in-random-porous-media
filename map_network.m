function network = map_network(domain, vertical_lines)
% to allow vertical lines: vertical_lines = 0, else vertical_lines = 1

image = domain.image;
X = domain.X*domain.pixel_scaling;
Y = domain.Y*domain.pixel_scaling;
r = domain.r*domain.pixel_scaling;

pores=[];
edgeX=[];
edgeY=[];
edge_length=[];
edge_width=[];

imagesizex=domain.size_x*domain.pixel_scaling;%size(image, 2);
imagesizey=domain.size_y*domain.pixel_scaling;%size(image, 1);
expanded_image=[image, image, image;
    image, image, image;
    image, image, image];


expanded_X=[X, X+imagesizex, X+2*imagesizex, X, X+imagesizex, X+2*imagesizex, X, X+imagesizex, X+2*imagesizex];
expanded_Y=[Y, Y, Y, Y+imagesizey, Y+imagesizey, Y+imagesizey, Y+2*imagesizey, Y+2*imagesizey, Y+2*imagesizey];
positions1=unique([expanded_X.', expanded_Y.'], 'rows'); %remove dublicate cylinder positions
expanded_X=positions1(:, 1).';
expanded_Y=positions1(:, 2).';


dt = delaunayTriangulation(expanded_X.', expanded_Y.');
[pore_sites, C]=voronoiDiagram(dt); %pore_sites = full list of all pore locations including those outside domain
n=length(C);


for i=1:n
    Vertices=[];
    Vert=pore_sites(C{i},:);%unique(, 'rows')
    for j=1:length(Vert)
        if isempty(Vertices)==1
            Vertices=[Vertices; Vert];
        else
            if ismember(Vert, Vertices, 'rows')==0
                Vertices=[Vertices; Vert];
            end
        end
    end
    x_centre=sum(Vertices(:, 1))/length(Vertices);
    y_centre=sum(Vertices(:, 2))/length(Vertices);
    if x_centre>imagesizex-2*r && x_centre<2*imagesizex+2*r && y_centre>imagesizey-2*r && y_centre<2*imagesizey+2*r
        pore_given=[Vertices(:, 1), Vertices(:, 2)];
        for j=1:size(pore_given, 1)
            if pore_given(j, 1)>=imagesizex && pore_given(j, 1)<=2*imagesizex && pore_given(j, 2)>=imagesizey && pore_given(j, 2)<=2*imagesizey
                if isempty(pores)==1
                    pores=[pore_given(j, :)];
                else
                    if ismember(pore_given(j, :), pores, 'rows')==0
                        pores=[pores; pore_given(j, :)];
                    end
                end
            end
        end
        
        %check if edges fall within [imagesize, 2*imagesize], else map into
        %correct area
        for j=1:size(pore_given, 1)
            if j==length(pore_given)
                posX=[pore_given(j, 1), pore_given(1, 1)];
                posY=[pore_given(j, 2), pore_given(1, 2)];
            else
                posX=[pore_given(j, 1), pore_given(j+1, 1)];
                posY=[pore_given(j, 2), pore_given(j+1, 2)];
            end
            if posX(1)>=imagesizex && posX(1)<=2*imagesizex && posX(2)>=imagesizex && posX(2)<=2*imagesizex && posY(1)>=imagesizey && posY(1)<=2*imagesizey && posY(2)>=imagesizey && posY(2)<=2*imagesizey  
                if isempty(edgeX)==1
                    k=0;
                    if posX(1)==posX(2)
                        if vertical_lines==1
                            k=k+1;
                        end
                        if round(posX(1), 3)==round(imagesizex, 3) || round(posX(1), 3)==round(2*imagesizex, 3)
                            k=k+1;
                        elseif posY(1)==posY(2)
                            k=k+1;
                        end
                    end
                    if posY(1)==imagesizey && posY(2)==imagesizey
                        k=k+1;
                    end
                    if k==0
                        edgeX=[edgeX; posX];
                        edgeY=[edgeY; posY];
                    end
                else
                    k=0;
                    if ismember(posX, edgeX, 'rows')==1
                        [~, indx]=ismember(posX, edgeX, 'rows');
                        if posY(1)==edgeY(indx, 1) && posY(2)==edgeY(indx, 2)
                            k=k+1;
                        end
                        if posX(1)==posX(2)
                            if vertical_lines == 1
                                k=k+1;
                            end
                            if ismember([posY(2), posY(1)], edgeY, 'rows')==1
                                k=k+1;
                            end
                            
                        end
                    elseif ismember([posX(2), posX(1)], edgeX, 'rows')==1
                        [~, indx]=ismember([posX(2), posX(1)], edgeX, 'rows');
                        if posY(2)==edgeY(indx, 1) && posY(1)==edgeY(indx, 2)
                            k=k+1;
                        end
                    end
                    if posX(1)==posX(2)
                        if vertical_lines == 1
                            k=k+1;
                        end
                        if round(posX(1), 3)==round(imagesizex, 3) || round(posX(1), 3)==round(2*imagesizex, 3)
                            k=k+1;
                        elseif posY(1)==posY(2)
                            k=k+1;
                        end
                    end
                    if posY(1)==imagesizey && posY(2)==imagesizey
                        k=k+1;
                    end
                    if k==0
                        edgeX=[edgeX; posX];
                        edgeY=[edgeY; posY];
                    end
                end
            else
                k=0;
                if posX(1)<=imagesizex || posX(1)>=2*imagesizex
                    m=(posY(2)-posY(1))/(posX(2)-posX(1));
                    c=posY(2)-m*posX(2);
                    if posX(2)<=imagesizex || posX(2)>=2*imagesizex
                        k=k+1;
                    end
                    if posX(1)<=imagesizex
                        posX(1)=imagesizex;
                    else
                        posX(1)=2*imagesizex;
                    end
                    posY(1)=m*posX(1)+c;
                elseif posX(2)<=imagesizex || posX(2)>=2*imagesizex
                    m=(posY(2)-posY(1))/(posX(2)-posX(1));
                    c=posY(1)-m*posX(1);
                    if posX(2)<=imagesizex
                        posX(2)=imagesizex;
                    else
                        posX(2)=2*imagesizex;
                    end
                    posY(2)=m*posX(2)+c;
                end
                if posX(1)==posX(2)
                    if vertical_lines ==1
                        k=k+1;
                    end
                    if round(posX(1), 3)==round(imagesizex, 3) || round(posX(1), 3)==round(2*imagesizex, 3)
                        k=k+1;
                    elseif posY(1)==posY(2)
                            k=k+1;
                    end
                end
                if posY(1)<=imagesizey || posY(1)>=2*imagesizey
                    if posX(1) == posX(2)
                        if vertical_lines==1
                            k=k+1;
                        end
                        if posY(2)<=imagesizey || posY(2)>=2*imagesizey
                            k=k+1;
                        end
                        if posY(1)<=imagesizey
                            posY(1)=imagesizey;
                        else
                            posY(1)=2*imagesizey;
                        end
                    else
                        m=(posY(2)-posY(1))/(posX(2)-posX(1));
                        c=posY(2)-m*posX(2);
                        if posY(2)<=imagesizey || posY(2)>=2*imagesizey
                            k=k+1;
                        end
                        if posY(1)<=imagesizey
                            posY(1)=imagesizey;
                        else
                            posY(1)=2*imagesizey;
                        end
                        posX(1)=(posY(1)-c)/m;
                    end
                elseif posY(2)<=imagesizey || posY(2)>=2*imagesizey
                    if posX(1) == posX(2)
                        if vertical_lines==1
                            k=k+1;
                        end
                        if posY(2)<=imagesizey
                            posY(2)=imagesizey;
                        else
                            posY(2)=2*imagesizey;
                        end
                    else
                        m=(posY(2)-posY(1))/(posX(2)-posX(1));
                        c=posY(1)-m*posX(1);
                        if posY(2)<=imagesizey
                            posY(2)=imagesizey;
                        else
                            posY(2)=2*imagesizey;
                        end
                        posX(2)=(posY(2)-c)/m;
                    end
                end
                if isempty(edgeX)==0
                    if ismember(posX, edgeX, 'rows')==1
                        if posX(1)==posX(2)
                            if vertical_lines==1
                                k=k+1;
                            end
                            indx1 = find(edgeX(:, 1)==posX(1));
                            indx2 = find(edgeX(:, 2)==posX(1));
                            indx = intersect(indx1, indx2);
                            if isempty(indx)==0
                                for ij=1:length(indx)
                                    if posY(1)==edgeY(indx(ij), 1) && posY(2)==edgeY(indx(ij), 2)
                                        k=k+1;
                                    elseif posY(2)==edgeY(indx(ij), 1) && posY(1)==edgeY(indx(ij), 2)
                                        k=k+1;
                                    end
                                end
                            end
                        else
                            [~, indx]=ismember(posX, edgeX, 'rows');
                            if posY(1)==edgeY(indx, 1) && posY(2)==edgeY(indx, 2)
                                k=k+1;
                            elseif posY(2)==edgeY(indx, 1) && posY(1)==edgeY(indx, 2)
                                k=k+1;
                            end
                        end
                    elseif ismember([posX(2), posX(1)], edgeX, 'rows')==1
                        [~, indx]=ismember([posX(2), posX(1)], edgeX, 'rows');
                        if posY(2)==edgeY(indx, 1) && posY(1)==edgeY(indx, 2)
                            k=k+1;
                        end
                    end
                end
                if posY(1)==imagesizey && posY(2)==imagesizey
                    k=k+1;
                end
                if k==0
                    edgeX=[edgeX; posX];
                    edgeY=[edgeY; posY];
                end
            end
        end
        
        %include pores on edge of domain in 'pores'
        for j=1:size(edgeX, 1)
            point1=[edgeX(j, 1), edgeY(j, 1)];
            point2=[edgeX(j, 2), edgeY(j, 2)];
            if isempty(pores)==1
                pores=[pores; point1];
            elseif ismember(point1, pores, 'rows')==0
                pores=[pores; point1];
            end
            if ismember(point2, pores, 'rows')==0
                pores=[pores; point2];
            end
        end
    end
end

problem_indices=[]; %check pores for NaN
for j=1:size(pores, 1)
    if isnan(pores(j, 1))==1 || isnan(pores(j, 2))==1
        problem_indices=[problem_indices, j];        
    end
end
for j=sort(problem_indices, 'descend')
    pores(j, :)=[];
end

problem_indices=[]; %check edges for NaN
for j=1:size(edgeX, 1)
    if isnan(edgeX(j, 1))==1 || isnan(edgeX(j, 2))==1 || isnan(edgeY(j, 1))==1 || isnan(edgeY(j, 2))==1
        problem_indices=[problem_indices, j];        
    end
end
for j=sort(problem_indices, 'descend')
    edgeX(j, :)=[];
    edgeY(j, :)=[];
end

pores=pores-[imagesizex*ones(length(pores), 1), imagesizey*ones(length(pores), 1)];
edgeX=edgeX-imagesizex;
edgeY=edgeY-imagesizey;

pores_interior1=[]; pores_interior2=[]; pores_interior3=[]; pores_interior4=[]; pores_exterior=[];
for i=1:length(pores) %re-order pores
    if pores(i, 1)==0
        pores_interior1=[pores_interior1; pores(i, :)];
    elseif pores(i, 2)==imagesizey
        pores_interior2=[pores_interior2; pores(i, :)];
    elseif pores(i, 1)==imagesizex
        pores_interior3=[pores_interior3; pores(i, :)];
    elseif pores(i, 2)==0
        pores_interior4=[pores_interior4; pores(i, :)];
    else
        pores_exterior=[pores_exterior; pores(i, :)];
    end
end
pores=[pores_interior1; pores_interior3; pores_interior2; pores_interior4; pores_exterior];

edgeX_holder1=[];
edgeX_holder2=[];
edgeY_holder1=[];
edgeY_holder2=[];
for i=1:size(edgeX, 1) %make it so edges go from left to right
    if edgeX(i, 1)>edgeX(i, 2)
        edgeX_holder1=[edgeX_holder1; edgeX(i, 2)];
        edgeX_holder2=[edgeX_holder2; edgeX(i, 1)];
        edgeY_holder1=[edgeY_holder1; edgeY(i, 2)];
        edgeY_holder2=[edgeY_holder2; edgeY(i, 1)];
    else
        edgeX_holder1=[edgeX_holder1; edgeX(i, 1)];
        edgeX_holder2=[edgeX_holder2; edgeX(i, 2)];
        edgeY_holder1=[edgeY_holder1; edgeY(i, 1)];
        edgeY_holder2=[edgeY_holder2; edgeY(i, 2)];
    end
end
edgeX=[edgeX_holder1, edgeX_holder2];
edgeY=[edgeY_holder1, edgeY_holder2];

edge_holder=unique([edgeX, edgeY], 'rows');
edgeX=edge_holder(:, 1:2);
edgeY=edge_holder(:, 3:4);

edgeX_i1=[]; edgeX_i2=[]; edgeX_i3=[]; edgeX_i4=[]; edgeX_e=[];
edgeY_i1=[]; edgeY_i2=[]; edgeY_i3=[]; edgeY_i4=[]; edgeY_e=[];

for i=1:size(edgeX, 1) %re-order edges
    if edgeX(i, 1)==0 || edgeX(i, 2)==0
        edgeX_i1=[edgeX_i1; edgeX(i, :)];
        edgeY_i1=[edgeY_i1; edgeY(i, :)];
    elseif edgeX(i, 1)==imagesizex || edgeX(i, 2)==imagesizex
        edgeX_i3=[edgeX_i3; edgeX(i, :)];
        edgeY_i3=[edgeY_i3; edgeY(i, :)];
    elseif edgeY(i, 1)==imagesizey || edgeY(i, 2)==imagesizey
        edgeX_i2=[edgeX_i2; edgeX(i, :)];
        edgeY_i2=[edgeY_i2; edgeY(i, :)];
    elseif edgeY(i, 1)==0 || edgeY(i, 2)==0
        edgeX_i4=[edgeX_i4; edgeX(i, :)];
        edgeY_i4=[edgeY_i4; edgeY(i, :)];
    else
        edgeX_e=[edgeX_e; edgeX(i, :)];
        edgeY_e=[edgeY_e; edgeY(i, :)];
    end
end

edgeX=[edgeX_i1; edgeX_i3; edgeX_i2;  edgeX_i4; edgeX_e];
edgeY=[edgeY_i1; edgeY_i3; edgeY_i2;  edgeY_i4; edgeY_e];

for j=1:length(edgeX) %calculate edge widths, respective radii and lengths
    if edgeY(j, 1)==edgeY(j, 2) && edgeY(j, 1)==0
        mex=(edgeX(j, 1)+edgeX(j, 2))/2;%find midpoint of edge
        mey=(edgeY(j, 1)+edgeY(j, 2))/2;
        [~, v1y]=closestpoint([mex, mey], [domain.X.', domain.Y.']);
        width = abs(v1y-r);
        edge_width=[edge_width; width];
    elseif edgeY(j, 1)==edgeY(j, 2) && edgeY(j, 1)==imagesizey
        mex=(edgeX(j, 1)+edgeX(j, 2))/2;%find midpoint of edge
        mey=(edgeY(j, 1)+edgeY(j, 2))/2;
        [~, v1y]=closestpoint([mex, mey], [domain.X.', domain.Y.']);
        width = abs(imagesizey - v1y-r);
        edge_width=[edge_width; width];
    else
        mex=(edgeX(j, 1)+edgeX(j, 2))/2+imagesizex;%find midpoint of edge
        mey=(edgeY(j, 1)+edgeY(j, 2))/2+imagesizey;
        [v1x, v1y]=closestpoint([mex, mey], [expanded_X.', expanded_Y.']);
        X_holder=expanded_X;
        Y_holder=expanded_Y;
        i1=find(expanded_X==v1x);
        if length(i1)>1
            A=i1;
            B=find(expanded_Y==v1y);
            i1=intersect(A, B);
        end
        X_holder(i1)=[];
        Y_holder(i1)=[];
        [v2x, v2y]=closestpoint([mex, mey], [X_holder.', Y_holder.']);
        il2=find(expanded_X==v2x);
        if length(il2)>1
            A=il2;
            B=find(expanded_Y==v2y);
            il2=intersect(A, B);
        end

        width=sqrt((v1x-v2x)^2+(v1y-v2y)^2)-2*r;
        edge_width=[edge_width; width];
    end
    len=sqrt((edgeX(j, 1)-edgeX(j, 2))^2+(edgeY(j, 1)-edgeY(j, 2))^2);
    edge_length=[edge_length; len];
end

Inc=zeros(size(edgeX, 1), length(pores));
for i=1:size(edgeX, 1) %determine incidence matrix
    for j=1:length(pores)
        if edgeX(i, 1)==pores(j, 1)
            if edgeY(i, 1)==pores(j, 2)
                Inc(i, j)=-1;
            elseif edgeX(i, 1) == edgeX(i, 2)
                if edgeY(i, 2)==pores(j, 2)
                    Inc(i, j)=1;
                end
            end
        elseif edgeX(i, 2)==pores(j, 1)
            if edgeY(i, 2)==pores(j, 2)
                Inc(i, j)=1;
            end
        end
    end
end

network.A = Inc;
network.pores = pores/domain.pixel_scaling;
network.edgeX = edgeX/domain.pixel_scaling;
network.edgeY = edgeY/domain.pixel_scaling;
network.edge_width = edge_width/domain.pixel_scaling;
network.edge_length = edge_length/domain.pixel_scaling;

network.M = length(edge_width);
network.K = length(pores);

k=0; k1=0; k2=0;
for i=1:network.K
    if pores(i, 1)==0
        k=k+1;
        k1=k1+1;
    elseif pores(i, 1)==imagesizex
        k=k+1;
        k2=k2+1;
    end
end
network.k = k; network.k1=k1; network.k2=k2;

m=0; m1=0; m2=0;
for i=1:network.M
    if edgeX(i, 1)==0
        m=m+1;
        m1=m1+1;
    elseif edgeX(i, 2)==imagesizex
        m=m+1;
        m2=m2+1;
    end
end
network.m = m; network.m1=m1; network.m2=m2;
end


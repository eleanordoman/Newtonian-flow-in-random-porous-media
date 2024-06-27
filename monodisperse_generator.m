function domain = monodisperse_generator(domain)
%Places circles of radius r until no more can be added where default is that no circles can touch edge of image

%min=minium distance between each circle
min = domain.min_width*domain.pixel_scaling;

%image=binary image of size imagesize_x x imagesize_y

imagesize_x = domain.size_x*domain.pixel_scaling;
imagesize_y  = domain.size_y*domain.pixel_scaling;

%X, Y=vector of centres of circles

r = domain.r*domain.pixel_scaling;

%min_boundary = minium distance between circles and edge of domain
min_boundary = domain.min_boundary*domain.pixel_scaling;
% take min_boundary<0 to allow circles to overlap with edges

tol=1000;

image = zeros(imagesize_y, imagesize_x);
[x, y] = meshgrid(1:imagesize_x, 1:imagesize_y);

if min_boundary>=0 %assume circles cannot be outside domain

    X=[randi([round(r+min_boundary), round(imagesize_x-r-min_boundary)], 1)]; %arrays to record centre of circles and first circle
    Y=[randi([round(r+min_boundary), round(imagesize_y-r-min_boundary)], 1)];
    counter=0;
    image((x-X(1)).^2+(y-Y(1)).^2<=r^2)=1; %creating first circle

    while counter<tol
        counter=0;
        overlap=1; %assume they overlap
        add=1;
        while (overlap>0.5)% && counter<tol %pick centres of circles until found set that do not overlap
            k=[];
            x1=randi([round(r+min_boundary), round(imagesize_x-(r+min_boundary))], 1);
            y1=randi([round(r+min_boundary), round(imagesize_y-(r+min_boundary))], 1);
            for j=1:length(X) %check that circle does not overlap other circles
                if (x1-X(j))^2+(y1-Y(j))^2<(2*r+min)^2
                    k=[k, 1];
                end 
            end 
        overlap=sum(k);
        counter=counter+1;
        if counter>tol
            add=0;
            break
        end
        end
        if add==1
            X=[X, x1]; %add circle centre coordinates to list of centres
            Y=[Y, y1];
            image((x-x1).^2+(y-y1).^2<=r^2)=1; %add circle to image
        end
    end
else %allow circles to be outside domain
    X=[randi([1, imagesize_x], 1)]; %arrays to record centre of circles and first circle
    Y=[randi([1, imagesize_y], 1)];
    counter=0;
    image((x-X(1)).^2+(y-Y(1)).^2<=r^2)=1; %creating first circle
    x1 = X; y1=Y;
    if x1<r
        if y1<r
            X = [X, x1+imagesize_x];
            Y = [Y, y1+imagesize_y];
            image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
        elseif y1>imagesize_y-r 
            X = [X, x1+imagesize_x];
            Y = [Y, y1-imagesize_y];
            image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
        else
            X = [X, x1+imagesize_x];
            Y = [Y, y1];
            image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
        end
    elseif x1>imagesize_x-r
        if y1<r
            X = [X, x1-imagesize_x];
            Y = [Y, y1+imagesize_y];
            image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
        elseif y1>imagesize_x-r
            X = [X, x1-imagesize_x];
            Y = [Y, y1-imagesize_y];
            image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
        else
            X = [X, x1-imagesize_x];
            Y = [Y, y1];
            image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
        end
    end
    if y1<r
        X = [X, x1];
        Y = [Y, y1+imagesize_y];
        image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
    elseif y1>imagesize_y-r
        X = [X, x1];
        Y = [Y, y1-imagesize_y];
        image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
    end
    while counter<tol
        counter=0;
        overlap=1; %assume they overlap
        add=1;
        while (overlap>0.5)% && counter<tol %pick centres of circles until found set that do not overlap
            k=[];
            x1=randi([1, imagesize_x], 1);
            y1=randi([1, imagesize_y], 1);
            for j=1:length(X) %check that circle does not overlap other circles
                if (x1-X(j))^2+(y1-Y(j))^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j)+imagesize_x)^2+(y1-Y(j))^2<(2*r+min)^2 %ensure no overlapping when cell is periodic
                    k=[k, 1];
                elseif (x1-X(j)-imagesize_x)^2+(y1-Y(j))^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j))^2+(y1-Y(j)+imagesize_y)^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j))^2+(y1-Y(j)-imagesize_y)^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j)+imagesize_x)^2+(y1-Y(j)-imagesize_y)^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j)+imagesize_x)^2+(y1-Y(j)+imagesize_y)^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j)-imagesize_x)^2+(y1-Y(j)-imagesize_y)^2<(2*r+min)^2
                    k=[k, 1];
                elseif (x1-X(j)-imagesize_x)^2+(y1-Y(j)+imagesize_y)^2<(2*r+min)^2
                    k=[k, 1];
                end 
            end 
        overlap=sum(k);
        counter=counter+1;
        if counter>tol
            add=0;
            break
        end
        end
        if add==1
            X=[X, x1]; %add circle centre coordinates to list of centres
            Y=[Y, y1];
            image((x-x1).^2+(y-y1).^2<=r^2)=1; %add circle to image
            x2 = x1; y2=y1;
            if x2<r
                X = [X, x2+imagesize_x];
                Y = [Y, y2];
                image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
            elseif x2>imagesize_x-r
                X = [X, x2-imagesize_x];
                Y = [Y, y2];
                image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
            end
            if y2<r
                X = [X, x2];
                Y = [Y, y2+imagesize_y];
                image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
            elseif y2>imagesize_y-r
                X = [X, x2];
                Y = [Y, y2-imagesize_y];
                image((x-X(end)).^2+(y-Y(end)).^2<=r^2)=1;
            end
        end
    end
end

domain.X=X/domain.pixel_scaling;
domain.Y=Y/domain.pixel_scaling;
domain.image=image;

end


function [X, Y] = closestpoint(P, A)
%find the closest point in A to P
if size(A, 1)==1
    X=A(1);
    Y=A(1);
else
    PX=P(1);
    PY=P(2);
    dist=zeros(length(A), 1);
    for i=1:length(dist)
        dist(i)=(A(i, 1)-PX)^2+(A(i, 2)-PY)^2;
    end
    [~, I]=min(dist);
    X=A(I, 1);
    Y=A(I, 2);
end
end


function newtonian_solve = calculate_vessel_resistance(domain, network, mu)
    
    edge_width = network.edge_width;
    depth = domain.depth;
    r = domain.r;
    edge_length = network.edge_length;
    
    if mu==0
        mu = 1.2*8.9*10^(-7);% default plasma value
    end
    
    n=length(edge_width);
    R = zeros(n, 1);
    for i=1:n
        if edge_width(i) == 0
            R(i)=0;
        elseif network.edgeY(i, 1)==network.edgeY(i, 2)
            if network.edgeY(i, 1)==0 || network.edgeY(i, 1)== domain.size_y
                R(i) = (sqrt(2)*mu*pi^4*r^0.5/(8*depth^2*edge_width(i)^1.5))*g(edge_width(i)/depth, edge_length(i)/(2*edge_width(i)*r)^0.5);
            else
                R(i) = (mu*pi^4*r^0.5/(8*depth^2*edge_width(i)^1.5))*g(edge_width(i)/depth, edge_length(i)/(edge_width(i)*r)^0.5);
            end
        else
            R(i) = (mu*pi^4*r^0.5/(8*depth^2*edge_width(i)^1.5))*g(edge_width(i)/depth, edge_length(i)/(edge_width(i)*r)^0.5);
        end
    end
    
    newtonian_solve.R = R;
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

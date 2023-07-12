function [X, S] = P_solver(D,rho,a)
[U, S, V] = svd(D,'econ');
[X, T1] = solve_P(S,rho,a,U,V);
S = diag(T1);
end
function [X,S1] = solve_P(S,rho,a,U,V)
    lambda = 1/rho;
    S0 = diag(S);
    m = length(S0);
    S1 = zeros(m,1);
    if lambda <= 1/(2*a^2)
        t = (lambda*a/2);
    else
        t = max(sqrt(2*lambda)-1/(2*a),0);
    end
    for i=1:m
        if S0(i) > t 
            S1(i) = prox(lambda,a,S0(i));
        else
            S1(i) = 0;
        end
    end
    X = U *diag(S1) *V';
end

function [s1] = prox(lambda,a,gamma)
    p1 = (1+a*abs(gamma))/3;
    f = phi(lambda,a,gamma);
    p2 = 1+2*cos(f/3 - pi/3);
    s1 = sign(gamma)*((p1*p2 - 1)/a);
end

function [f] = phi(lambda,a,gamma)
    num = 27 * lambda * (a^2);
    den = 2 * (1 + a*abs(gamma))^3;
    f = acos(num/den - 1);
end
   


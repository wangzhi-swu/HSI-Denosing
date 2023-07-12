function [L,S] = NFF_ALM(X,opts,a,C)
[m,n] = size(X);

if isfield(opts,'lambda')
    lambda = opts.lambda;
else
    lambda = 6e-3;
end

if isfield(opts,'type')
    type = opts.type;
else
    type = 1;
end

if isfield(opts,'mu')
   mu = opts.mu;
else
    mu = 8e-3;
end

if isfield(opts,'rate')
   rate = opts.rate;
else
    rate = 1.5;
end

if isfield(opts,'tol')
   tol = opts.tol;
else
    tol = 1e-3;
end

if isfield(opts,'max_itr')
    max_itr = opts.max_itr;
else
    max_itr = 500;
end
S = zeros(m,n);  % sparse noise
Lambda = S;      % Lagrangian multiplier
% mu =  min(2/(0.99*norm(X)+(1/(2*a)))^2,a/(0.99*norm(X)));
mu = 2/(norm(X))^2;
% fprintf('mu = %f \n',mu);
norm_X = norm(X,'fro');
for iter = 1:max_itr
    temp = X - S - Lambda/mu; % 
    % step 1: Solving P(L)
    [L,ss] = P_solver( temp,mu,a);
     % step 2: Solving S
    [S] = S_solver(Lambda,X,L,lambda,mu,type); % 
     % step 3: Solving Lambda
    Lambda = Lambda + mu*(L + S - X); % 
    Rel_Err = norm (X - L - S ,'fro')/norm_X; % convergence condition
    mu = mu* rate; % rho 
    if Rel_Err < tol        
        break;
    end   
end
% disp(['  relative error',num2str(Rel_Err)]);
end



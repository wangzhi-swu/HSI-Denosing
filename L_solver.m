function [ X, T] = L_solver(D,rho,T0,a,strname)
% D = Xk - Sk - Lambda/mu
% rho = mu/2
% T0 = sig
% a = gamma
% strname
% solve argmin_L λ||L||_{2,1} + (ρk/2)||S-Wk||_F^2
[U, S, V] = svd(D,'econ'); % svd of D

for t = 1:100
   [X, T1] = DCInner(S,rho,T0,a,U,V,strname);
   err = sum((T1-T0).^2);
   if err < 1e-6
       break
   end
   T0 = T1;
end
T = T1;
end

function [X,t] = DCInner(S,rho,J,epislon,U,V,strname)
   % S = Σ
   % rho = rho
   % J = T0 = sig
   % epislon = a = gamma
   % U V strname
lambda = 0.5/rho;
S0 = diag(S); % Σii
   if strcmp(strname,'Lap')
       grad = exp (-J/epislon)/epislon;% ▼Φ(sig(i)) % Laplace 
    else 
       grad = (1+epislon)*epislon./(epislon+J).^2; % Geman
    end

t = max(S0 - lambda * grad,0); % WSVT ： max(Σii - (▼Φ(sig(i)))/ρk，0)
X = U *diag(t) * V'; % L*
end


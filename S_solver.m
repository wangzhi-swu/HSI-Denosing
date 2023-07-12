function [error]= S_solver(Y1,X,W,lambda,mu,type)
% Y1 = Lambda
% X
% W = L
% mu = rho
% lambda type
% solve argmin_S = λ||S||_{2,1} + (ρk/2)||S-Wk||_F^2
switch type
    case 1 % l1范数
        D = -Y1/mu+(X-W);
        E = zeros(size(D));
        epsilon = lambda/mu;
        DD = abs(D)-epsilon;
        DD2 = DD.*sign(D);
        ID = abs(D)>epsilon;
        E(ID) = DD2(ID);
    case 21 % l2,1范数
        alpha=lambda/mu; 
        G=X-W-Y1/mu; % G : Wk in Lemma2
        G1 = sqrt(sum(G.^2,1)); % ||Wk(:,j)||_2
        G1(G1==0) = alpha;
        G2 = (G1-alpha)./G1;
        E = G*diag((G1>alpha).*G2);
end
error = E;
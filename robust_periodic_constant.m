function [K,niu,rho, Qr, Sr, Rr] = robust_periodic_constant(X,U,R,n,L,N,niu,rho)


%% Solving data driven problem

cvx_solver mosek

cvx_begin sdp

%  cvx_solver_settings('maxit',150)
 variable P(n,n,N+1) symmetric
 variable Y(L,n,N)
 variable Qr(n,n) symmetric
 variable Sr(n,L) 
 variable Rr(L,L) symmetric

%  variable niu 
%  variable rho
 % niu = 1;
 % rho = 450;
%  niu = 1;
%  rho = 450; %CHANGE BACK TO 230
%  

Rr-eps*eye(L) <= 0;
%Rr = -gamm*eye(L);
[eye(n);R']'*[Qr Sr;Sr' Rr]*[eye(n);R']>=0;

for i = 1:N

    [P(:,:,i+1)-eye(n)-Qr,-Sr,X(:,:,i+1)*Y(:,:,i);
        -Sr',-Rr,Y(:,:,i);
        Y(:,:,i)'*X(:,:,i+1)',Y(:,:,i)',P(:,:,i)] -eps*eye(2*n+L)>= 0;

    X(:,:,i)*Y(:,:,i) == P(:,:,i);
    P(:,:,i) - niu*eye(n) >= 0;
    P(:,:,i) - rho*eye(n) <= 0;


end
P(:,:,1) == P(:,:,N+1);
cvx_end
 
 
%% Calculating the gain 

for i = 1 : N
    
    K(:,:,i) = U(:,:,i)*Y(:,:,i)*pinv(P(:,:,i));
end

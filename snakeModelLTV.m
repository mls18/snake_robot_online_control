function [dx_ltv_fun_sim, A_ltv]=snakeModelLTV(ct,cn,ct_exp,cn_exp,N)


%% Snake parameters values

A = zeros(N-1, N);
D = zeros(N-1, N);

for i = 1: N-1
    
    A(i,i)= 1;
    A(i,i+1) = 1;
    
    D(i,i) = 1;
    D(i, i+1) = -1;
end

D_bar = D'*inv(D*D'); 
e_bar = ones(N-1, 1); 
m = 1; %kg
l = 0.14; %m
J = 0.0016; %kgm^2
lamda_1 = 0.5;
lamda_2 = 20;
cp=(cn-ct)/(2*l);
%% Reference trajectory parameters
alpha = 0.045;
omega = 120*pi/180;
delta = 40*pi/180;
offset = 0;

%% Reference trajectory
syms t a omg d gam
for j = 1:N-1
    
    phi(j,:) = a*sin(omg*t +(j-1)*d) + gam;
    phi_dot(j,:) = omg*a*cos(omg*t +(j-1)*d);
    phi_dotdot(j,:)=-omg^2*a*sin(omg*t +(j-1)*d);
    
end
phi_ref = subs(phi,[a, omg, d, gam],[alpha, omega, delta, offset]);
phi_dot_ref = subs(phi_dot,[a, omg, d, gam],[alpha, omega, delta, offset]);
phi_dotdot_ref = subs(phi_dotdot,[a, omg, d, gam],[alpha, omega, delta, offset]);


phi_ref_fun = matlabFunction(phi_ref);
dphi_ref_fun = matlabFunction(phi_dot_ref);
ddphi_ref_fun = matlabFunction(phi_dotdot_ref);

%% Errors in coefficient of frictions (=0 when exact)

cp_exp = (cn_exp-ct_exp)/(2*l);  % previoulsy around 5 cp + 0.05


%% LTV error coordinate model
syms t 
d_x = sym('d_x',[2*N,1]);
d_u = sym('d_u',[N-1,1]);

dx(1:2*N-2) = [zeros(N-1) eye(N-1);(cp-cp_exp)/m*d_x(2*N-1)*A*D' (cn_exp-cn)/m*eye(N-1)]*d_x(1:2*N-2) + [zeros(N-1); eye(N-1)]*d_u;
dx(2*N-1:2*N) = [-ct/m, 2*cp/(N*m)*e_bar'*phi_ref_fun(t);2*cp/(N*m)*e_bar'*phi_ref_fun(t), -cn/m]*d_x(2*N-1:2*N)+[-cp/(N*m)*phi_ref_fun(t)'*A*D_bar*dphi_ref_fun(t);0];
dx_ltv_fun_sim = matlabFunction(dx', 'vars', {t, d_x, d_u});
A_ltv = [zeros(N-1) eye(N-1);(cp-cp_exp)/m*d_x(2*N-1)*A*D' (cn_exp-cn)/m*eye(N-1)];
A_ltv = matlabFunction(A_ltv);
end
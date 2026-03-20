function [dx_fun_sim,phi_ref_fun,dphi_ref_fun,ddphi_ref_fun,d_sys_fun,cn_exp,cp_exp] = snakeModel_noisy(ct,cn,ct_exp,cn_exp,N,noise_switch)
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

     
%% Snake model function
x = sym('x',[2*(2*N+6),1]);
u_sys = sym('u_sys',[N-1,1]);

num_original_states = length(x)/2;  % Assuming we add one noise state per original state

noise_states = x(num_original_states+1:end);

% Calculate derivative of noise states (first-order low-pass filter)
tau_a = 0.2;  % Time constant of the filter
sigma_a = 0.002;       % Noise amplitude

tau_v = 0.05;
sigma_v = 0.01;

tau_p = 0.1;
sigma_p = 0.001;

tau = [tau_a*ones(N,1);tau_p*ones(2,1);tau_a*ones(N-1,1);tau_v*ones(5,1)];
sigma = [sigma_a*ones(N,1);sigma_p*ones(2,1);sigma_a*ones(N-1,1);sigma_v*ones(5,1)]
% Each noise state follows: dx_noise = -x_noise/tau + sigma*sqrt(2/tau)*randn
dx_noise = noise_switch*(-noise_states./tau + (sigma.*sqrt(2./tau)).*randn(size(noise_states)));
   
dx(1:N-1) = x(N+3:2*N+1);
dx(N) = x(2*N+2);
dx(N+1) = x(2*N+3)*cos(x(N)) - x(2*N+4)*sin(x(N));
dx(N+2) = x(2*N+3)*sin(x(N)) + x(2*N+4)*cos(x(N));
dx(N+3:2*N+1) = u_sys + ddphi_ref_fun(t) -cn/m*x(N+3:2*N+1)+cp/m*x(2*N+3)*A*D'*x(1:N-1)+cn_exp/m*(x(N+3:2*N+1)+x(3*N+9:4*N+7))-cp_exp/m*(x(2*N+3)+x(4*N+9))*A*D'*(x(1:N-1)+x(2*N+7:3*N+5));
dx(2*N+2) = -lamda_1*x(2*N+2) + lamda_2/(N-1)*x(2*N+3)*e_bar'*x(1:N-1);
dx(2*N+3) = -ct/m*x(2*N+3) +2*cp/(N*m)*x(2*N+4)*e_bar'*x(1:N-1) -cp/(N*m)*x(1:N-1)'*A*D_bar*x(N+3:2*N+1);
dx(2*N+4) = -cn/m*x(2*N+4) + 2*cp/(N*m)*x(2*N+3)*e_bar'*x(1:N-1);
dx(2*N+5) = -ct/m*x(2*N+5) +2*cp/(N*m)*x(2*N+6)*e_bar'*phi_ref_fun(t) -cp/(N*m)*phi_ref_fun(t)'*A*D_bar*dphi_ref_fun(t);
dx(2*N+6) = -cn/m*x(2*N+6) + 2*cp/(N*m)*x(2*N+5)*e_bar'*phi_ref_fun(t);
dx = [dx';dx_noise];

dx_fun_sim = matlabFunction(dx, 'vars', {t, x, u_sys});

d_sys = [zeros(N-1,1);(cp-cp_exp)/m*(x(2*N+3)-x(2*N+5))*A*D'*(x(1:N-1)-phi_ref_fun(t))] + [zeros(N-1,1); (cn_exp-cn)/m*(dphi_ref_fun(t)) + (cp-cp_exp)/m*x(2*N+5)*A*D'*(phi_ref_fun(t))];
d_sys_fun = matlabFunction(d_sys);

end
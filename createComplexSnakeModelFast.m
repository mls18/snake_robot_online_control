function [dx_fun_sim, phi_ref_func, dphi_ref_func, ddphi_ref_func, n] = ...
    createComplexSnakeModelFast(c_t, c_n, c_t_tilde, c_n_tilde, N)
%% Creates complex snake robot model with specified friction coefficients
%
% KEY CHANGE vs. original: fully numeric — no Symbolic Toolbox, no matlabFunction.
% Constant matrices are precomputed once. The ODE RHS evaluates numerically
% at each step, using backslash (\) for the Schur complement solve instead
% of symbolic inv(). Scales to N=7+ without issue.
%
% Inputs:
%   c_t, c_n           : actual tangential and normal friction coefficients
%   c_t_tilde, c_n_tilde: nominal friction coefficients (for PFL)
%   N                  : number of links
% Outputs:
%   dx_fun_sim         : @(t, x, u_sys) dynamics handle
%   phi_ref_func, dphi_ref_func, ddphi_ref_func: reference trajectory handles
%   n                  : number of links (compatibility)

n = N;

%% Physical parameters
m_link = 1;
l      = 0.14;
J_link = (m_link * l^2) / 3;

%% ---- Constant matrices (all numeric) ------------------------------------

% Addition matrix A  (n-1)×n
A = zeros(n-1, n);
for i = 1:n-1
    A(i, i)   = 1;
    A(i, i+1) = 1;
end

% Difference matrix D  (n-1)×n
D = zeros(n-1, n);
for i = 1:n-1
    D(i, i)   =  1;
    D(i, i+1) = -1;
end

e    = ones(n, 1);
O    = (1:2:2*n-1)';          % [1, 3, 5, ..., 2n-1]

DDt_inv  = inv(D * D');        % (n-1)×(n-1), purely numeric
N_mat    = A' * DDt_inv * D;   % n×n
V_mat    = A' * DDt_inv * A;   % n×n

% Mapping R: theta = R * [phi_1 ... phi_{n-1}, theta_N]'
R_mat = zeros(n, n);
for i = 1:n
    R_mat(i, :) = [zeros(1, i-1), ones(1, n-i+1)];
end

l_N      = l * N_mat;
l_N_tr   = l_N';
m_lsq_V  = m_link * l^2 * V_mat;

%% ---- Reference trajectory (analytic, no syms) ---------------------------
alpha_ref  = 13.9 * pi / 180;
omega_ref  = 120  * pi / 180;
delta_ref  = 40  * pi / 180;
offset_ref = 0;
jvec = (0:n-2)';              % (n-1)×1 phase-offset vector

phi_ref_func   = @(t)  alpha_ref * sin(omega_ref * t + jvec * delta_ref) + offset_ref;
dphi_ref_func  = @(t)  omega_ref * alpha_ref * cos(omega_ref * t + jvec * delta_ref);
ddphi_ref_func = @(t) -omega_ref^2 * alpha_ref * sin(omega_ref * t + jvec * delta_ref);

%% ---- Pack parameters for the numeric dynamics closure -------------------
p.n          = n;
p.m          = m_link;
p.J          = J_link;
p.c_t        = c_t;
p.c_n        = c_n;
p.c_t_tilde  = c_t_tilde;
p.c_n_tilde  = c_n_tilde;
p.R          = R_mat;
p.e          = e;
p.l_N_tr     = l_N_tr;
p.l_N        = l_N;
p.m_lsq_V    = m_lsq_V;
p.ddphi_ref  = ddphi_ref_func;

dx_fun_sim = @(t, x, u_sys) complexSnakeDyn(t, x, u_sys, p);

end % createComplexSnakeModel


%% ==========================================================================
function dx = complexSnakeDyn(t, x, u_sys, p)
%% Numeric evaluation of complex snake dynamics at each ODE step.
%
% State layout (length 2*(n+2)):
%   x(1   : n-1 )  = phi_1 ... phi_{n-1}
%   x(n   : n+2 )  = [theta_N, p_x, p_y]
%   x(n+3 : 2n+1)  = phi_dot
%   x(2n+2: 2n+4)  = [theta_N_dot, p_x_dot, p_y_dot]

n = p.n;

%% Absolute joint angles and velocities
theta     = p.R * x(1:n);
theta_dot = p.R * x(n+3:2*n+2);

S = diag(sin(theta));   % n×n
C = diag(cos(theta));   % n×n

%% Link CM velocities in global frame
x_dot_links = p.l_N_tr * (S * theta_dot) + p.e * x(2*n+3);
y_dot_links = -p.l_N_tr * (C * theta_dot) + p.e * x(2*n+4);
z_dot = [x_dot_links; y_dot_links];   % 2n×1

%% Friction forces (actual and nominal)
f       = computeFriction(p.c_t,       p.c_n,       S, C, z_dot);
f_tilde = computeFriction(p.c_t_tilde, p.c_n_tilde, S, C, z_dot);

fx  = f(1:n);        fy  = f(n+1:2*n);
fxt = f_tilde(1:n);  fyt = f_tilde(n+1:2*n);

%% Inertia matrix M and Coriolis matrix W
M_n = p.J * eye(n) + S * p.m_lsq_V * S + C * p.m_lsq_V * C;
W_n = S * p.m_lsq_V * C - C * p.m_lsq_V * S;

%% Augmented inertia M_bar  (n+2)×(n+2)
RtMR  = p.R' * M_n * p.R;
M_bar = [RtMR,        zeros(n, 2); ...
         zeros(2, n), p.n * p.m * eye(2)];

% Partition: block-1 = joint DOFs (1:n-1), block-2 = [theta_N, px, py] (n:n+2)
M11 = M_bar(1:n-1, 1:n-1);   % (n-1)×(n-1)
M12 = M_bar(1:n-1, n:end);   % (n-1)×3
M21 = M_bar(n:end, 1:n-1);   % 3×(n-1)
M22 = M_bar(n:end, n:end);   % 3×3  (cheap numeric inverse)

M22_inv = inv(M22);

%% Coriolis: W*diag(v)*v == W*(v.*v)
W_bar = [p.R' * W_n * (theta_dot .* theta_dot); zeros(2, 1)];
W2    = W_bar(n:end);   % 3×1

%% G matrices (joint torques and CM forces from friction)
S_lN = S * p.l_N;
C_lN = C * p.l_N;
G_bar = [-p.R' * S_lN,  p.R' * C_lN; ...
         -p.e',          zeros(1, n); ...
          zeros(1, n),  -p.e'        ];
G1 = G_bar(1:n-1, :);   % (n-1)×2n
G2 = G_bar(n:end, :);   % 3×2n

%% PFL components
K_pfl = -M22_inv * (W2 + G2 * [fx; fy]);   % 3×1
Jm    = -M22_inv * M21;                      % 3×(n-1)

%% Friction error (nominal vs. actual)
df         = [fxt - fx; fyt - fy];                      % 2n×1
term_error = (G1 - M12 * M22_inv * G2) * df;            % (n-1)×1

%% *** Key fix: numeric backslash — no symbolic inv of (n-1)×(n-1) matrix ***
M_schur = M11 - M12 * M22_inv * M21;                    % (n-1)×(n-1)
accel   = u_sys + p.ddphi_ref(t) + M_schur \ term_error; % (n-1)×1

%% State derivative
dx = [x(n+3:2*n+1);         % phi_dot          (n-1)×1
      x(2*n+2:2*n+4);       % [θ_N, px, py] dot  3×1
      accel;                 % phi_ddot         (n-1)×1
      K_pfl + Jm * accel];  % [θ_N, px, py] ddot 3×1

end % complexSnakeDyn


%% ==========================================================================
function f = computeFriction(c_t, c_n, S, C, z_dot)
%% Viscous friction force vector for given coefficients and link velocities
F = -[c_t * C*C + c_n * S*S,    (c_t - c_n) * S*C; ...
      (c_t - c_n) * S*C,         c_n * C*C + c_t * S*S];
f = F * z_dot;
end % computeFriction
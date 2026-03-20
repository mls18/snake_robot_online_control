%% main_adaptive_complex.m
% Identical to main_adaptive2 but uses the complex model for data collection.
% The simplified LTV (snakeModelLTV) is kept for D[k] and A_d[k].
% Because the complex model has no embedded noise states, V[k] = 0
% and R[k] = A_d[k]*0 - 0 - D[k] = -D[k].

clear; clc

%% -----------------------------------------------------------------------
%% Parameters  (unchanged)
%% -----------------------------------------------------------------------
N         = 5;
n         = (N-1)*2;
m         = N-1;
omega     = 120*pi/180;
alpha     = 2*pi/omega;
Ts        = 0.005;
N_period  = round(alpha/Ts);
N_sim     = 20000;

rho_sdp       = 8000;
eta_sdp       = 1;
safety_factor = 5;

L             = n + m;
k_d           = 5;
k_i           = 20;
sigma_explore = 6.0;

A_mat = zeros(N-1, N);
D_mat = zeros(N-1, N);
for i = 1:N-1
    A_mat(i,i) = 1;  A_mat(i,i+1) = 1;
    D_mat(i,i) = 1;  D_mat(i,i+1) = -1;
end
m_kg = 1;
l    = 0.14;

%% -----------------------------------------------------------------------
%% Initial friction coefficients  (unchanged)
%% -----------------------------------------------------------------------
cn_nom  = 3.2;  ct_nom  = 1.2;
cn_true = 3;    ct_true = 1;

%% -----------------------------------------------------------------------
%% PART 1: Offline controller design
%% -----------------------------------------------------------------------
fprintf('=================================================================\n');
fprintf('PART 1: Offline controller design\n');
fprintf('=================================================================\n\n');

%% --- CHANGE: use complex model instead of snakeModel_noisy ---
[dx_fun_sim, phi_ref_fun, dphi_ref_fun, ddphi_ref_fun] = ...
    createComplexSnakeModelFast(ct_true, cn_true, ct_nom, cn_nom, N);

%% LTV unchanged - still snakeModelLTV for D[k] and A_d[k]
[dx_ltv_fun_sim, ~] = snakeModelLTV(ct_true, cn_true, ct_nom, cn_nom, N);

[A_d, vx_ref, vy_ref] = computeAd(dx_ltv_fun_sim, cn_nom, ct_nom, ...
                                    cn_true, ct_true, N, n, m_kg, l, ...
                                    Ts, N_period, A_mat, D_mat);

%% Offline data collection
fprintf('Collecting offline data (%d experiments x %d steps)...\n', L, N_period);
N_steps   = N_period * L;
U_explore = randn(m, L, N_period) * sigma_explore;

X_out = zeros(n, N_steps+1);
Z_out = zeros(n, N_steps+1);
U_out = zeros(m, N_steps);
D_out = zeros(n, N_steps);
V_out = 0.001*randn(n,N_steps+1);

%% --- CHANGE: complex model state is 2*N+4 (no noise state block) ---
x0        = 0.01 * rand(2*N+4, 1);


X_out(:,1) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(0); dphi_ref_fun(0)];
Z_out(:,1) = X_out(:,1);

for i = 1:N_steps
    k_per = mod(i-1, N_period) + 1;
    j_exp = ceil(i / N_period);

    u_i        = k_d*(dphi_ref_fun((i-1)*Ts) - x0(N+3:2*N+1)) ...
               + k_i*(phi_ref_fun((i-1)*Ts)  - x0(1:N-1)) ...
               + U_explore(:, j_exp, k_per);
    U_out(:,i) = u_i;

    %% --- CHANGE: integrate complex model ---
    [~, x_nl] = ode45(@(t,x) dx_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0);
    x0 = x_nl(end,:)';

    X_out(:,i+1) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)];
    Z_out(:,i+1) = X_out(:,i+1)+ V_out(:,i+1);

    % D[k]: one-step LTV prediction error (LTV reinitialised from true state)
    x0_ltv = [X_out(:,i); vx_ref(k_per); vy_ref(k_per)];
    [~, x_lv] = ode45(@(t,x) dx_ltv_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0_ltv);
    D_out(:,i) = X_out(:,i+1) - x_lv(end, 1:n)';
end

%% Reshape into 3-D data matrices
Z = zeros(n, L, N_period+1);
U = zeros(m, L, N_period);
V = zeros(n, L, N_period+1);   % zeros throughout
D = zeros(n, L, N_period);
for j = 1:L
    idx      = N_period*(j-1)+1 : N_period*j;
    Z(:,j,:) = Z_out(:, [idx, idx(end)+1]);
    U(:,j,:) = U_out(:, idx);
    D(:,j,:) = D_out(:, idx);
    V(:,j,:) = V_out(:, [idx, idx(end)+1]);
    
end

%% R[k] = A_d[k]*V[k] - V[k+1] - D[k]  =>  R[k] = -D[k]  (V=0)
R = zeros(n, L, N_period);
for k = 1:N_period
    R(:,:,k) = A_d(:,:,k)*V(:,:,k) - V(:,:,k+1) - D(:,:,k);
end

%% Calculation max R
sigma_max_R = zeros(1, N);
for k = 1:N
    sigma_max_R(k) = norm(R(:,:,k), 2);
end
[rho_R,k_idx] = max(sigma_max_R);


% Construct R_bar with sigma_max = rho_R via SVD definition
U_rand = orth(randn(n, n));        % random orthogonal matrix in R^{n x n}
V_rand = orth(randn(L, L));        % random orthogonal matrix in R^{L x L}

Sigma = zeros(n, L);
Sigma(1,1) = rho_R;                % only one non-zero singular value = rho_R

R_bar = U_rand * Sigma * V_rand';  % sigma_max(R_bar) = rho_R exactly

fprintf('Solving offline SDP...\n');
tic
[K,niu,rho, Qr, Sr, Rr] = robust_periodic_constant(Z,U,R_bar,n,L,N_period,eta_sdp,rho_sdp);

fprintf('Offline SDP solved in %.1f s\n\n', toc);

%% -----------------------------------------------------------------------
%% PART 2: Online closed-loop simulation
%% -----------------------------------------------------------------------
fprintf('=================================================================\n');
fprintf('PART 2: Online simulation\n');
fprintf('=================================================================\n\n');

%% --- CHANGE: complex model state size ---
x0        = rand(2*N+4, 1);


x_error      = zeros(N_sim+1, n);
t_sampled    = zeros(1, N_sim+1);
vt           = zeros(1, N_sim+1);
x_error(1,:) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(0); dphi_ref_fun(0)];
%% --- CHANGE: vt has no noise state offset in complex model ---
vt(1)        = x0(2*N+3);
t = [];  x = [];

a = 30;
epsi1 = 1e-7;
epsi2 = 5e-6;
mu     = N_period*2;

gc      = 0;
g_count = 0;
k_gc    = 0;
k_exp   = Inf;
j_col   = 1;

N_collect = N_period * L;
X_col = zeros(n, N_collect+1);
Z_col = zeros(n, N_collect+1);
V_col = 0.001*randn(n,N_sim+1);  
U_col = zeros(m, N_collect);
D_col = zeros(n, N_collect);
U_explore_online = [];

for i = 1:N_sim

    i_c = mod(i-1, N_period) + 1;

    % Ground change
    if i == 3000
        fprintf('[sim] Ground change at step %d (t = %.2f s)\n', i, (i-1)*Ts);
        cn_true = 16;  ct_true = 4;
        %% --- CHANGE: rebuild complex model ---
        [dx_fun_sim, ~, ~, ~] = createComplexSnakeModelFast(ct_true, cn_true, ct_nom, cn_nom, N);
    end

    % Control input
    zeta = x0([1:N-1, N+3:2*N+1])+ V_col(:,i) - [phi_ref_fun((i-1)*Ts); dphi_ref_fun((i-1)*Ts)];
    if  i >= k_exp && i <= k_exp + N_collect - 1 && ~isempty(U_explore_online)
        k_per_col = mod(j_col-1, N_period) + 1;
        j_exp_col = min(ceil(j_col / N_period), L);
        u_i = K(:,:,i_c) * zeta + U_explore_online(:, j_exp_col, k_per_col);
    else
        u_i = K(:,:,i_c) * zeta;
    end

    %% --- CHANGE: integrate complex model ---
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
    [t_out, x_out] = ode45(@(t,x) dx_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0,opts);
    x0 = x_out(end,:)';

    if i < N_sim
        t = [t; t_out(1:end-1)];
        x = [x; x_out(1:end-1,:)];
    else
        t = [t; t_out];
        x = [x; x_out];
    end
    x_error(i+1,:) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)];
    %% --- CHANGE: vt index ---
    vt(i+1) = (cos(x_out(end, N+2)) * x_out(end, 2*N+3) + sin(x_out(end, N+2)) * x_out(end, 2*N+4));
    

    t_sampled(i+1) = t_out(end);

    % Ground change detection
    if i > 2*a + 6*pi/(omega*Ts) + 1 && gc == 0 && g_count == 0
        s = detectGround(vt, a, epsi1, epsi2, pi/(omega*Ts), i+1);
        if s == 1
            fprintf('[detect] Ground change detected at step %d (t = %.2f s)\n', i+1, i*Ts);
            k_gc = i + 1;
            gc   = 1;
        end
    end
    
    % Nominal update
     if gc == 1 && i == k_gc + mu + 1
         fprintf('[update] Updating nominal at step %d (t = %.2f s)...\n', i, (i-1)*Ts);

         u_bar    = ddphi_ref_fun((i-1)*Ts) + u_i;
         x_a(1,:) = x_error(i-1,:) + [phi_ref_fun((i-1)*Ts); dphi_ref_fun((i-1)*Ts)]'+ V_col(:,i-1)';
         x_a(2,:) = x_error(i,:)   + [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)]'+ V_col(:,i)';
         x_a(3,:) = x_error(i+1,:) + [phi_ref_fun((i+1)*Ts); dphi_ref_fun((i+1)*Ts)]'+ V_col(:,i+1)';  

         [cn_nom_u, ct_nom_u] = updateNominal(x_a, u_bar, vt(k_gc:k_gc+mu)*0.6253/mean(vt(k_gc:k_gc+mu)), ...
             ct_nom, cn_nom, 1, l, Ts, N);

        
       
         if abs(cn_nom_u - cn_true) <= 2
             cn_nom = cn_nom_u;
         else cn_nom = 16.5;
         end

         if abs(ct_nom_u - ct_true) <= 2
            ct_nom = ct_nom_u;
         else ct_nom = 4.5;
         end
         
        %% --- CHANGE: rebuild complex model with updated nominal ---
        [dx_fun_sim, ~, ~, ~] = createComplexSnakeModelFast(ct_true, cn_true, ct_nom, cn_nom, N);
        [dx_ltv_fun_sim, ~]   = snakeModelLTV(ct_true, cn_true, ct_nom, cn_nom, N);

        [A_d, vx_ref, vy_ref] = computeAd(dx_ltv_fun_sim, cn_nom, ct_nom, ...
                                            cn_true, ct_true, N, n, m_kg, l, ...
                                            Ts, N_period, A_mat, D_mat);

        U_explore_online = randn(m, L, N_period)*sigma_explore;

        k_exp        = i + 1;
        j_col        = 1;
        X_col(:,1)   = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)];
        Z_col(:,1)   = X_col(:,1)+V_col(:,1);
        g_count      = 1;
    end

    % Online data collection
    if gc == 1 && i >= k_exp && i <= k_exp + N_collect - 1
        k_per = mod(j_col-1, N_period) + 1;

        U_col(:,j_col)   = u_i;
        X_col(:,j_col+1) = x_error(i+1,:)';
        Z_col(:,j_col+1) = X_col(:,j_col+1)+ V_col(:,j_col+1);

        x0_ltv     = [X_col(:,j_col); vx_ref(k_per); vy_ref(k_per)];
        [~, x_lv]  = ode45(@(t,x) dx_ltv_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0_ltv);
        D_col(:,j_col) = X_col(:,j_col+1) - x_lv(end, 1:n)';

        j_col = j_col + 1;
    end

    % Controller redesign
    if gc == 1 && i == k_exp + N_collect - 1
        fprintf('[redesign] Data collection complete. Solving online SDP...\n');

        Z_new = zeros(n, L, N_period+1);
        U_new = zeros(m, L, N_period);
        V_new = zeros(n, L, N_period+1);   % zeros
        D_new = zeros(n, L, N_period);
        for jj = 1:L
            idx           = N_period*(jj-1)+1 : N_period*jj;
            Z_new(:,jj,:) = Z_col(:, [idx, idx(end)+1]);
            U_new(:,jj,:) = U_col(:, idx);
            D_new(:,jj,:) = D_col(:, idx);
            V_new(:,jj,:) = V_col(:, [idx, idx(end)+1]);
        end

        R_new = zeros(n, L, N_period);
        for k = 1:N_period
            R_new(:,:,k) = A_d(:,:,k)*V_new(:,:,k) - V_new(:,:,k+1) - D_new(:,:,k);
        end
        sigma_max_R = zeros(1, N);
        for k = 1:N
            sigma_max_R(k) = norm(R_new(:,:,k), 2);
        end
        [rho_R,k_idx] = max(sigma_max_R);

        % Construct R_bar with sigma_max = rho_R via SVD definition
        U_rand = orth(randn(n, n));        % random orthogonal matrix in R^{n x n}
        V_rand = orth(randn(L, L));        % random orthogonal matrix in R^{L x L}

        Sigma = zeros(n, L);
        Sigma(1,1) = rho_R;                % only one non-zero singular value = rho_R
        
        R_bar_new = 0.1*U_rand * Sigma * V_rand';
        %R_bar_new = 0.1*R_new(:,:,k_idx);  % sigma_max(R_bar) = rho_R exactly

    
        tic

        [K,niu,rho, Qr_new, Sr_new, Rr_new] = robust_periodic_constant(Z_new,U_new,R_bar_new,n,L,N_period,eta_sdp,rho_sdp);
        fprintf('[redesign] Online SDP solved in %.1f s\n\n', toc);

        gc      = 0;

        j_col   = 1;
        k_exp   = Inf;
    end

end

%% -----------------------------------------------------------------------
%% Plots  (unchanged)
%% -----------------------------------------------------------------------
figure; hold on; grid on
stairs(t_sampled, x_error(:,1:N-1), 'LineWidth', 1.2)
xline(3000*Ts, 'r--', 'LineWidth', 1.5)
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('$\delta\phi$ [m]', 'Interpreter', 'latex', 'FontSize', 13)
title('Joint angle tracking error', 'FontSize', 12)
legend(arrayfun(@(i) sprintf('$\\delta\\phi_%d$', i), 1:N-1, 'UniformOutput', false), ...
       'Interpreter', 'latex')

figure; hold on; grid on
z1_norm = vecnorm(x_error, 2, 2);
k_steps = 0:N_sim;
decay   = sqrt(rho/niu) * (1 - 1/rho).^(k_steps/2) * z1_norm(1);
plot(t_sampled, z1_norm, 'b',   'LineWidth', 1.5, 'DisplayName', '$\|z_1[k]\|$')
plot(t_sampled, decay,   'g--', 'LineWidth', 1.5, 'DisplayName', '$\sqrt{\rho/\eta}(1-1/\rho)^{k/2}\|z_1[0]\|$')
xline(3000*Ts, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ground change')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('$\|z_1[k]\|$', 'Interpreter', 'latex', 'FontSize', 13)
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 11)
title('Actuated error norm vs decaying bound', 'FontSize', 12)

% Actual joint angles for joints 1 and 3 vs reference
%phi_ref_sampled = zeros(N_sim+1, N-1);
for k = 1:N_sim+1
    phi_ref_sampled(k,:) = phi_ref_fun(t_sampled(k))';
end
phi_actual = x_error(:,1:N-1) + phi_ref_sampled;

figure; hold on; grid on
plot(t_sampled, phi_actual(:,1),      'b',   'LineWidth', 1.2, 'DisplayName', '$\phi_1$ actual')
plot(t_sampled, phi_ref_sampled(:,1), 'b--', 'LineWidth', 1.0, 'DisplayName', '$\phi_1$ reference')
% plot(t_sampled, phi_actual(:,3),      'r',   'LineWidth', 1.2, 'DisplayName', '$\phi_3$ actual')
% plot(t_sampled, phi_ref_sampled(:,3), 'r--', 'LineWidth', 1.0, 'DisplayName', '$\phi_3$ reference')
xline(3000*Ts, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ground change')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Joint angle [rad]', 'Interpreter', 'latex', 'FontSize', 13)
title('Actual vs reference joint angles ($\phi_1$, $\phi_3$)', 'Interpreter', 'latex', 'FontSize', 12)
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 10)
grid on
%% -----------------------------------------------------------------------
%% Local function: computeAd  (unchanged)
%% -----------------------------------------------------------------------
function [A_d, vx_ref, vy_ref] = computeAd(dx_ltv_fun_sim, cn_nom, ct_nom, ...
                                             cn_true, ct_true, N, n, m_kg, l, ...
                                             Ts, N_period, A_mat, D_mat)
    u_zero        = zeros(N-1, 1);
    x0_ref        = zeros(2*N, 1);
    x0_ref(2*N-1) = 0.0048;
    x0_ref(2*N)   = -8.058461988272277e-05;

    vx_ref    = zeros(1, N_period+1);
    vy_ref    = zeros(1, N_period+1);
    vx_ref(1) = x0_ref(2*N-1);
    vy_ref(1) = x0_ref(2*N);

    for k = 1:N_period
        [~, xr]     = ode45(@(t,x) dx_ltv_fun_sim(t, x, u_zero), ...
                             [(k-1)*Ts, k*Ts], x0_ref);
        x0_ref      = xr(end,:)';
        vx_ref(k+1) = x0_ref(2*N-1);
        vy_ref(k+1) = x0_ref(2*N);
    end

    cp_nom  = (cn_nom  - ct_nom)  / (2*l);
    cp_true = (cn_true - ct_true) / (2*l);
    A_d     = zeros(n, n, N_period);
    for k = 1:N_period
        B21        = (cp_nom - cp_true)/m_kg * vx_ref(k) * A_mat * D_mat';
        Ac         = [zeros(N-1), eye(N-1); ...
                      B21, (cn_true - cn_nom)/m_kg * eye(N-1)];
        A_d(:,:,k) = eye(n) + Ac * Ts;
    end
end
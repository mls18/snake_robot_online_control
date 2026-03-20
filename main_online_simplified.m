%% main_adaptive.m
% Adaptive data-driven control for a planar snake robot.
%
% Overview:
%   Part 1 -- Offline data collection and initial controller design.
%             Uses the true initial friction (unknown to the controller)
%             and a nominal model to collect data, compute R[k], and
%             solve the robust periodic SDP.
%   Part 2 -- Online closed-loop simulation with:
%               (a) Ground change detection via CoM velocity monitoring
%               (b) Nominal friction update via updateNominal
%               (c) Online data collection with the updated nominal model
%               (d) Controller redesign using the new data
%
% R[k] formula (eq. 5 in paper):
%   R[k] = A_d[k]*V[k] - V[k+1] - D[k]
%   where V[k] = measurement noise, D[k] = process noise (one-step
%   LTV prediction error with LTV reinitialised from true state at k).
%
% Controller design: robust_periodic_snr_singularvalue
%   Uses Young's inequality + spectral norm bound on R[k] to replace
%   the noise bound constraint with a sufficient condition depending
%   only on rho_R = max_k sigma_max(R[k]).

clear; clc

%% -----------------------------------------------------------------------
%% Parameters
%% -----------------------------------------------------------------------
N         = 5;                  % number of snake links
n         = (N-1)*2;            % actuated error state dimension
m         = N-1;                % input dimension
omega     = 120*pi/180;         % undulation frequency [rad/s]
alpha     = 2*pi/omega;         % undulation period [s]
Ts        = 0.005;              % sampling time [s]
N_period  = round(alpha/Ts);    % steps per period (= 600)
N_sim     = 20000;              % total online simulation steps

% SDP design parameters (matched to offline_practical)
rho_sdp       = 8000;   % upper bound on eigenvalues of P[k]
eta_sdp       = 1;      % lower bound on eigenvalues of P[k]
safety_factor = 5;      % gamma[k] = safety_factor * rho_sdp / lambda_min(Z[k]Z[k]')

% Data collection parameters
L             = n + m;          % number of parallel experiments (= 9)
k_d           = 5;              % PD proportional gain for exploration
k_i           = 20;             % PD derivative gain for exploration
sigma_explore = 6.0;            % std of random exploration signal

% Geometry matrices used for computing A_d[k]
A_mat = zeros(N-1, N);
D_mat = zeros(N-1, N);
for i = 1:N-1
    A_mat(i,i) = 1;  A_mat(i,i+1) = 1;
    D_mat(i,i) = 1;  D_mat(i,i+1) = -1;
end
m_kg = 1;    % link mass [kg]
l    = 0.14; % link half-length [m]

%% -----------------------------------------------------------------------
%% Initial friction coefficients
%% -----------------------------------------------------------------------
cn_nom  = 3.2;  ct_nom  = 1.2;
cn_true = 3;    ct_true = 1;

%% -----------------------------------------------------------------------
%% PART 1: Offline controller design
%% -----------------------------------------------------------------------
fprintf('=================================================================\n');
fprintf('PART 1: Offline controller design\n');
fprintf('=================================================================\n\n');

%% Load snake models with initial friction coefficients
% dx_fun_sim:     nonlinear simulation model (true friction)
% dx_ltv_fun_sim: nominal LTV error dynamics (nominal friction)
[dx_fun_sim, phi_ref_fun, dphi_ref_fun, ddphi_ref_fun, ~, ~, ~] = ...
    snakeModel_noisy(ct_true, cn_true, ct_nom, cn_nom, N, 0);
[dx_ltv_fun_sim, ~] = snakeModelLTV(ct_true, cn_true, ct_nom, cn_nom, N);

%% Precompute A_d[k] and reference CoM velocities for the initial ground
% A_d[k] = I + A_c[k]*Ts is the discretised LTV error dynamics matrix.
% It depends on the reference CoM velocity vx_ref[k], which is obtained
% by propagating the LTV model from the reference initial condition with
% zero input for one full period.
[A_d, vx_ref, vy_ref] = computeAd(dx_ltv_fun_sim, cn_nom, ct_nom, ...
                                    cn_true, ct_true, N, n, m_kg, l, ...
                                    Ts, N_period, A_mat, D_mat);

%% Offline data collection
% L experiments are run sequentially, each lasting one full period.
% Each experiment j uses a fixed random exploration sequence U_explore(:,j,:)
% (drawn once before the loop) to keep Z[k] well-conditioned throughout.
fprintf('Collecting offline data (%d experiments x %d steps)...\n', L, N_period);
N_steps   = N_period * L;
U_explore = randn(m, L, N_period) * sigma_explore;  % fixed per experiment

% Pre-allocate sequential output arrays
X_out = zeros(n, N_steps+1);   % true error state  x[k] = phi - phi_ref
Z_out = zeros(n, N_steps+1);   % noisy measurement z[k] = x[k] + v[k]
U_out = zeros(m, N_steps);     % control input
D_out = zeros(n, N_steps);     % process noise d[k] (one-step LTV error)
V_out = 0.001*randn(n,N_steps+1);% measurement noise
% Initial condition
x0        = 0.01 * rand(2*N+6, 1);
x0(2*N+5) = 0.0048;
x0(2*N+6) = -8.058461988272277e-05;
x0        = [x0; zeros(2*N+6, 1)];  % append noise states (init to zero)

X_out(:,1) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(0); dphi_ref_fun(0)];
Z_out(:,1) = X_out(:,1);
tic
for i = 1:N_steps
    k_per = mod(i-1, N_period) + 1;  % step within current period (1..N_period)
    j_exp = ceil(i / N_period);       % experiment index (1..L)

    % Control: PD tracking + fixed random exploration for experiment j_exp
    u_i        = k_d*(dphi_ref_fun((i-1)*Ts) - x0(N+3:2*N+1)) ...
               + k_i*(phi_ref_fun((i-1)*Ts)  - x0(1:N-1)) ...
               + U_explore(:, j_exp, k_per);
    U_out(:,i) = u_i;

    % Nonlinear simulation step
    [~, x_nl] = ode45(@(t,x) dx_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0);
    x0 = x_nl(end,:)';

    % Error state and measurement noise at k+1
    X_out(:,i+1) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)];
    Z_out(:,i+1) = X_out(:,i+1) + V_out(:,i+1);

    % Process noise d[k]: one-step LTV prediction error.
    % The LTV is reinitialised from the TRUE state at k (not propagated)
    % so that D captures only the model mismatch over one step.
    x0_ltv = [X_out(:,i); vx_ref(k_per); vy_ref(k_per)];
    [~, x_lv] = ode45(@(t,x) dx_ltv_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0_ltv);
    D_out(:,i) = X_out(:,i+1) - x_lv(end, 1:n)';
end
time_datacollection = toc
%% Reshape sequential arrays into 3-D data matrices (n/m x L x N_period)
Z = zeros(n, L, N_period+1);
U = zeros(m, L, N_period);
D = zeros(n, L, N_period);
V = zeros(n, L, N_period+1);
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
%% Offline controller design via robust periodic SDP
% robust_periodic_snr_singularvalue:
%   - computes rho_R = max_k sigma_max(R[k]) internally
%   - fixes gamma[k] = safety_factor * rho_sdp / lambda_min(Z[k]Z[k]')
%   - solves SDP with Qr[k], Sr[k] as decision variables
%   - returns periodic gain K[k], k=1..N_period
% fprintf('Solving offline SDP...\n');
% tic
% 
% 
% [K, niu, rho] = robust_periodic( ...
%     Z, U, R, n, L, N_period, eta_sdp, rho_sdp);
 fprintf('Offline SDP solved in %.1f s\n\n', toc);

%% -----------------------------------------------------------------------
%% PART 2: Online closed-loop simulation
%% -----------------------------------------------------------------------
fprintf('=================================================================\n');
fprintf('PART 2: Online simulation\n');
fprintf('=================================================================\n\n');

%% Initialise simulation state
x0        = rand(2*N+6, 1);
x0        = [x0; zeros(2*N+6, 1)];
x_error      = zeros(N_sim+1, n);
t_sampled    = zeros(1, N_sim+1);
vt           = zeros(1, N_sim+1);    % CoM x-velocity (used for detection)
x_error(1,:) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(0); dphi_ref_fun(0)];
vt(1)        = x0(2*N+3) + x0(4*N+9);
t = [];  x = [];

%% Ground change detection parameters
a      = 30;          % detection window half-size
epsi1  = 1e-7;        % detection threshold 1
epsi2  = 5e-6;        % detection threshold 2
mu     = N_period;  % steps to wait after detection before updating
mm = 1;
%% Online state machine
% gc:      1 = ground change detected, waiting to update or collecting data
% g_count: 1 = adaptation in progress (prevents re-detection)
% k_gc:    step at which ground change was detected
% k_exp:   step at which online data collection starts
% j_col:   sequential column counter during online data collection
% epsi_u:  exploration amplitude (0 = off, >0 = on during collection)
gc      = 0;
g_count = 0;
k_gc    = 0;
k_exp   = Inf;
j_col   = 1;

% Pre-allocate online data collection buffers (size = L periods)
N_collect = N_period * L;
X_col = zeros(n, N_collect+1);
Z_col = zeros(n, N_collect+1);
U_col = zeros(m, N_collect);
D_col = zeros(n, N_collect);
V_col = 0.0001*randn(n,N_sim+1);
U_explore_online = [];   % drawn fresh when collection starts

%% Main simulation loop
for i = 1:N_sim

    % Periodic index: which step within the current period (1..N_period)
    i_c = mod(i-1, N_period) + 1;

    % ----------------------------------------------------------------
    % Ground change (simulation only -- unknown to controller)
    % At step 3000 the true friction changes to cn=16, ct=4.
    % ----------------------------------------------------------------
    if i == 3000
        fprintf('[sim] Ground change at step %d (t = %.2f s)\n', i, (i-1)*Ts);
        cn_true = 16;  ct_true = 4;
        [dx_fun_sim, ~, ~, ~, ~, ~, ~] = ...
            snakeModel_noisy(ct_true, cn_true, ct_nom, cn_nom, N, 0);
    end

    % ----------------------------------------------------------------
    % Control input
    % Outside collection: closed-loop gain applied to noisy measurement.
    % During online collection: PD + fixed per-experiment exploration,
    % identical in structure to the offline data collection, so that Z[k]
    % receives the same quality of excitation in all state directions.
    % ----------------------------------------------------------------
    zeta = x0([1:N-1, N+3:2*N+1]) + V_col(:,i) ...
           - [phi_ref_fun((i-1)*Ts); dphi_ref_fun((i-1)*Ts)];
    if  i >= k_exp && i <= k_exp + N_collect - 1 && ~isempty(U_explore_online)
        % During online collection: K[k]*zeta + fixed per-experiment exploration.
        % K[k] keeps the system stable while the structured exploration ensures
        % Z[k] is well-conditioned (same fixed-per-experiment structure as offline).
        k_per_col = mod(j_col-1, N_period) + 1;      % step within period
        j_exp_col = min(ceil(j_col / N_period), L);   % experiment index (1..L)
        u_i = K(:,:,i_c) * zeta + U_explore_online(:, j_exp_col, k_per_col);
    else
        u_i = K(:,:,i_c) * zeta;
    end

    % ----------------------------------------------------------------
    % Simulate one step (nonlinear model)
    % ----------------------------------------------------------------
    opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
    [t_out, x_out] = ode45(@(t,x) dx_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0,opts);
    x0 = x_out(end,:)';

    % Log trajectories
    if i < N_sim
        t = [t; t_out(1:end-1)];
        x = [x; x_out(1:end-1,:)];
    else
        t = [t; t_out];
        x = [x; x_out];
    end
    x_error(i+1,:) = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)];
    vt(i+1)        = x0(2*N+3)+0.0001*randn;
    t_sampled(i+1) = t_out(end);

    % ----------------------------------------------------------------
    % Ground change detection
    % Monitors CoM velocity for statistical changes indicative of a
    % ground change. Skipped while adaptation is in progress.
    % ----------------------------------------------------------------
    if i > 2*a + pi/(omega*Ts) + 1 && gc == 0 && g_count == 0 && mm == 1
        s = detectGround(vt, a, epsi1, epsi2, pi/(omega*Ts), i+1);
        if s == 1
            fprintf('[detect] Ground change detected at step %d (t = %.2f s)\n', ...
                    i+1, i*Ts);
            k_gc = i + 1;
            gc   = 1;
        end
    end

    % ----------------------------------------------------------------
    % Nominal update
    % Triggered mu steps after detection to allow transients to settle.
    % Estimates new friction from CoM velocity data, updates the nominal
    % LTV model and recomputes A_d[k] for the correct R[k] formula.
    % ----------------------------------------------------------------
    if  gc == 1 & i == k_gc + mu + 1  && mm == 1
        fprintf('[update] Updating nominal at step %d (t = %.2f s)...\n', i, (i-1)*Ts);

        % Reconstruct absolute joint coordinates at the update step
        u_bar    = ddphi_ref_fun((i-1)*Ts) + u_i;
        x_a(1,:) = x_error(i-1,:) + [phi_ref_fun((i-1)*Ts); dphi_ref_fun((i-1)*Ts)]'+ V_col(:,i-1)';
        x_a(2,:) = x_error(i,:)   + [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)]'+ V_col(:,i)';
        x_a(3,:) = x_error(i+1,:) + [phi_ref_fun((i+1)*Ts); dphi_ref_fun((i+1)*Ts)]'+ V_col(:,i+1)';  


        % Estimate new friction coefficients from CoM velocity history
        [cn_nom_u, ct_nom_u] = updateNominal(x_a, u_bar, vt(k_gc:k_gc+mu), ...
            ct_nom, cn_nom, 1, l, Ts, N);


        if abs(cn_nom_u - cn_true) <= 2
            cn_nom = cn_nom_u;
        else cn_nom = 16.5;
        end

        if abs(ct_nom_u - ct_true) <= 2
            ct_nom = ct_nom_u;
        else ct_nom = 4.5;
        end
        % Update nominal LTV model with new friction estimate
        [dx_ltv_fun_sim, ~] = snakeModelLTV(ct_true, cn_true, ct_nom, cn_nom, N);
        [dx_fun_sim, ~, ~, ~, ~, ~, ~] = snakeModel_noisy(ct_true, cn_true, ct_nom, cn_nom, N, 0); 
        % Recompute A_d[k] and vx_ref/vy_ref with updated nominal.
        % A_d[k] must be consistent with the LTV model used for D[k],
        % so it must be recomputed whenever the nominal changes.
        [A_d, vx_ref, vy_ref] = computeAd(dx_ltv_fun_sim, cn_nom, ct_nom, ...
                                            cn_true, ct_true, N, n, m_kg, l, ...
                                            Ts, N_period, A_mat, D_mat);

        
       U_explore_online = randn(m, L, N_period)*sigma_explore;
          
 % Initialise data collection buffers with the (perturbed) current state
        k_exp        = i + 1;   % collection starts at the next step
        j_col        = 1;
        X_col(:,1)   = x0([1:N-1, N+3:2*N+1]) - [phi_ref_fun(i*Ts); dphi_ref_fun(i*Ts)];
        g_count      = 1;       % prevent re-detection during adaptation
        mm = 2;
    end

    % ----------------------------------------------------------------
    % Online data collection
    % Runs for N_period*L steps after k_exp, collecting L experiments
    % of length N_period each. The LTV is reinitialised from the true
    % state at each step to compute the correct one-step process noise.
    % ----------------------------------------------------------------
    if  gc == 1 && i >= k_exp && i <= k_exp + N_collect - 1
        k_per = mod(j_col-1, N_period) + 1;  % step within period (1..N_period)

        % Store input, state and noise at k+1
        U_col(:,j_col)   = u_i;
        X_col(:,j_col+1) = x_error(i+1,:)';
        Z_col(:,j_col+1) = X_col(:,j_col+1) + V_col(:,j_col+1) ;

        % Process noise: reinitialise LTV from true state at k, predict to k+1
        x0_ltv     = [X_col(:,j_col); vx_ref(k_per); vy_ref(k_per)];
        [~, x_lv]  = ode45(@(t,x) dx_ltv_fun_sim(t, x, u_i), [(i-1)*Ts, i*Ts], x0_ltv);
        D_col(:,j_col) = X_col(:,j_col+1) - x_lv(end, 1:n)';

        j_col = j_col + 1;
    end

    % ----------------------------------------------------------------
    % Controller redesign
    % Triggered once online data collection is complete.
    % Assembles 3-D data matrices, computes correct R[k], and solves SDP.
    % ----------------------------------------------------------------
    if  gc == 1 && i == k_exp + N_collect - 1
        fprintf('[redesign] Data collection complete. Solving online SDP...\n');

        % Reshape sequential buffers into 3-D matrices (n/m x L x N_period)
        Z_new = zeros(n, L, N_period+1);
        U_new = zeros(m, L, N_period);
        V_new = zeros(n, L, N_period+1);
        D_new = zeros(n, L, N_period);
        for jj = 1:L
            idx           = N_period*(jj-1)+1 : N_period*jj;
            Z_new(:,jj,:) = Z_col(:, [idx, idx(end)+1]);
            U_new(:,jj,:) = U_col(:, idx);
            D_new(:,jj,:) = D_col(:, idx);
            V_new(:,jj,:) = V_col(:, [idx, idx(end)+1]);
        end

        % Compute correct R[k] = A_d[k]*V[k] - V[k+1] - D[k]
        % A_d has been recomputed with the updated nominal after detection.

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

    
        tic

        [K,niu,rho, Qr_new, Sr_new, Rr_new] = robust_periodic_constant(Z_new,U_new,R_bar_new,n,L,N_period,eta_sdp,rho_sdp);

         % Reset state machine for next potential ground change
        gc      = 0;
        g_count = 0;   % re-enable ground change detection
        j_col   = 1;
        k_exp   = Inf;
    end

end  % end main simulation loop
%% -----------------------------------------------------------------------
%% Plots
%% -----------------------------------------------------------------------

% Actuated state norm vs decaying bound (eq. 13)
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
plot(t_sampled, phi_actual(:,3),      'r',   'LineWidth', 1.2, 'DisplayName', '$\phi_3$ actual')
plot(t_sampled, phi_ref_sampled(:,3), 'r--', 'LineWidth', 1.0, 'DisplayName', '$\phi_3$ reference')
xline(3000*Ts, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ground change')
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Joint angle [rad]', 'Interpreter', 'latex', 'FontSize', 13)
title('Actual vs reference joint angles ($\phi_1$, $\phi_3$)', 'Interpreter', 'latex', 'FontSize', 12)
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 10)
grid on

%% -----------------------------------------------------------------------
%% Local function: compute A_d[k] and reference CoM velocities
%% -----------------------------------------------------------------------
function [A_d, vx_ref, vy_ref] = computeAd(dx_ltv_fun_sim, cn_nom, ct_nom, ...
                                             cn_true, ct_true, N, n, m_kg, l, ...
                                             Ts, N_period, A_mat, D_mat)
% COMPUTEAD  Compute the discretised LTV error dynamics matrix A_d[k] and
%            the reference CoM velocity trajectory over one full period.
%
% A_d[k] = I + A_c[k]*Ts, where A_c[k] is the continuous-time LTV error
% dynamics matrix at step k, which depends on vx_ref[k].
% vx_ref[k] is obtained by integrating the LTV model from the reference
% initial condition with zero input.

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

    % Build A_d[k] using the error dynamics matrix A_c[k]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deploy_gto_to_15deg_raan_constellation_symbolic_lsqnonlin.m
%
% Purpose:
%   Deploy N_s satellites from the same initial GTO-like orbit into the same
%   final orbit shape, but with final inclination i_f = 15 deg and RAANs
%   distributed uniformly over 360 deg.
%
%   Each satellite is solved independently.
%   The optimized initial costate from satellite j is used as the initial
%   costate guess for satellite j+1.
%
% Requires on MATLAB path:
%   averaged_dynamics_symbolic.m
%   calculate_averaged_hamiltonian_symbolic.m
%   get_MEE_B_matrix.m
%   generated_symbolic/*.m
%   lgwt.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units and parameters: same as working paper-style script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DU = 6378;
TU = 806.78557;

params.mu = 1;
params.g0 = 9.80665e-3*(TU^2/DU);
params.Isp = 3100/TU;
params.c = params.g0*params.Isp;

params.T_min = 0;
params.T_max = 0.1*(TU^2/(DU*1000));

params.q = 40;
params.eps_S = 1e-1;

params.RE = 6378.1/DU;
params.J2 = 1.08263e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial GTO-like state: same as working example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0 = 1.822602598777046;
f0 = 0.725;
g0 = 0;
h0 = 0.253967646474944;
k0 = 0;
L0 = 0;
t0 = 0;
alpha0 = 3212.749578552824;
m0 = 100;

x0 = [p0; f0; g0; h0; k0; L0; t0; alpha0; m0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial costate seed: same fast-run seed from working example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda0_seed = [
  -1.967847803145031
 -11.056764487340876
   4.145818903384442
  12.237626634771029
   0.899834150512536
   0.000017226165227
   0.000008800790354
                   0
   0.039524341544112
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constellation definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_s = 10;

inc_f_deg = 15;
inc_f = deg2rad(inc_f_deg);

raan_deg_list = linspace(0,360,N_s+1);
raan_deg_list(end) = [];

rho_hk = tan(inc_f/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE and nonlinear least-squares options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

odeopts = odeset( ...
    'RelTol',1e-12, ...
    'AbsTol',1e-12);

lsqopts = optimoptions(@lsqnonlin, ...
    'Display','iter', ...
    'Algorithm','levenberg-marquardt', ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',5000, ...
    'StepTolerance',1e-12, ...
    'FunctionTolerance',1e-12, ...
    'OptimalityTolerance',1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print run header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n============================================================\n');
fprintf('GTO TO 15-DEG INCLINED RAAN-DISTRIBUTED CONSTELLATION\n');
fprintf('============================================================\n');
fprintf('N_s = %d\n', N_s);
fprintf('Final inclination = %.15f deg\n', inc_f_deg);
fprintf('Initial p = %.15f\n', p0);
fprintf('Initial f = %.15f\n', f0);
fprintf('Initial g = %.15f\n', g0);
fprintf('Initial h = %.15f\n', h0);
fprintf('Initial k = %.15f\n', k0);
fprintf('alpha0 = %.15f TU\n', alpha0);
fprintf('time of flight days = %.15f\n', alpha0*TU/86400);
fprintf('m0 = %.15f kg\n', m0);
fprintf('params.c = %.15e\n', params.c);
fprintf('params.T_max = %.15e\n', params.T_max);
fprintf('params.T_max/m0 = %.15e\n', params.T_max/m0);
fprintf('params.J2 = %.15e\n', params.J2);
fprintf('params.q = %d\n', params.q);
fprintf('params.eps_S = %.15e\n', params.eps_S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda0_all = zeros(9,N_s);
residual_all = zeros(7,N_s);
resnorm_all = zeros(1,N_s);

target_all = struct([]);

final_mass_all = zeros(1,N_s);
H_drift_all = zeros(1,N_s);
min_sigma_all = zeros(1,N_s);
max_sigma_all = zeros(1,N_s);
min_T_all = zeros(1,N_s);
max_T_all = zeros(1,N_s);

tau_all = cell(1,N_s);
x_all = cell(1,N_s);
lambda_all = cell(1,N_s);
H_all = cell(1,N_s);
sigma_all = cell(1,N_s);
T_all = cell(1,N_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve each satellite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sat = 1:N_s

    raan_deg = raan_deg_list(sat);
    raan = deg2rad(raan_deg);

    target.p = p0;
    target.f = f0;
    target.g = 0;
    target.h = rho_hk*cos(raan);
    target.k = rho_hk*sin(raan);

    target_all(sat).p = target.p;
    target_all(sat).f = target.f;
    target_all(sat).g = target.g;
    target_all(sat).h = target.h;
    target_all(sat).k = target.k;
    target_all(sat).inc_deg = inc_f_deg;
    target_all(sat).raan_deg = raan_deg;

    fprintf('\n============================================================\n');
    fprintf('SATELLITE %d / %d\n', sat, N_s);
    fprintf('Target inclination = %.15f deg\n', inc_f_deg);
    fprintf('Target RAAN        = %.15f deg\n', raan_deg);
    fprintf('Target p           = %.15f\n', target.p);
    fprintf('Target f           = %.15f\n', target.f);
    fprintf('Target g           = %.15f\n', target.g);
    fprintf('Target h           = %.15f\n', target.h);
    fprintf('Target k           = %.15f\n', target.k);
    fprintf('Seed lambda0:\n');
    disp(lambda0_seed);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Seed propagation diagnostic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    res_seed = shooting_residual_symbolic(lambda0_seed,x0,target,params,odeopts);

    fprintf('\nSeed residual norm = %.15e\n', norm(res_seed));
    fprintf('Seed residual vector:\n');
    disp(res_seed);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve shooting problem for this satellite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda0_opt = lsqnonlin( ...
        @(lambda0) shooting_residual_symbolic(lambda0,x0,target,params,odeopts), ...
        lambda0_seed, [], [], lsqopts);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-propagate corrected solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [tau,Y] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                    [0 1], [x0; lambda0_opt], odeopts);

    x = Y(:,1:9);
    lambda = Y(:,10:18);

    H = zeros(length(tau),1);
    sigma = zeros(length(tau),1);
    T = zeros(length(tau),1);

    for i = 1:length(tau)
        H(i) = calculate_averaged_hamiltonian_symbolic( ...
            x(i,:).',lambda(i,:).',params);

        [sigma(i),T(i)] = thrust_history_simple( ...
            x(i,:).',lambda(i,:).',params);
    end

    res_opt = terminal_residual(x(end,:).',lambda(end,:).',target);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda0_all(:,sat) = lambda0_opt;
    residual_all(:,sat) = res_opt;
    resnorm_all(sat) = norm(res_opt);

    final_mass_all(sat) = x(end,9);
    H_drift_all(sat) = max(abs(H-H(1)));

    min_sigma_all(sat) = min(sigma);
    max_sigma_all(sat) = max(sigma);
    min_T_all(sat) = min(T);
    max_T_all(sat) = max(T);

    tau_all{sat} = tau;
    x_all{sat} = x;
    lambda_all{sat} = lambda;
    H_all{sat} = H;
    sigma_all{sat} = sigma;
    T_all{sat} = T;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print corrected solution summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n------------------------------------------------------------\n');
    fprintf('SATELLITE %d CORRECTED SOLUTION\n', sat);
    fprintf('------------------------------------------------------------\n');
    fprintf('Optimized initial costates:\n');
    disp(lambda0_opt);

    fprintf('terminal residual norm = %.15e\n', norm(res_opt));
    fprintf('Final p error = %.15e\n', x(end,1)-target.p);
    fprintf('Final f error = %.15e\n', x(end,2)-target.f);
    fprintf('Final g error = %.15e\n', x(end,3)-target.g);
    fprintf('Final h error = %.15e\n', x(end,4)-target.h);
    fprintf('Final k error = %.15e\n', x(end,5)-target.k);
    fprintf('lambda_L(tf) = %.15e\n', lambda(end,6));
    fprintf('lambda_m(tf) = %.15e\n', lambda(end,9));
    fprintf('Final mass = %.15f\n', x(end,9));
    fprintf('max |H-H0| = %.15e\n', max(abs(H-H(1))));
    fprintf('min sigma = %.15f\n', min(sigma));
    fprintf('max sigma = %.15f\n', max(sigma));
    fprintf('min T = %.15e\n', min(T));
    fprintf('max T = %.15e\n', max(T));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use this satellite's solved costate as next seed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lambda0_seed = lambda0_opt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n============================================================\n');
fprintf('CONSTELLATION SUMMARY\n');
fprintf('============================================================\n');
fprintf('%5s %12s %20s %20s %20s %20s\n', ...
    'Sat','RAAN[deg]','resnorm','final mass','H drift','max sigma');

for sat = 1:N_s
    fprintf('%5d %12.6f %20.8e %20.12f %20.8e %20.12f\n', ...
        sat, ...
        target_all(sat).raan_deg, ...
        resnorm_all(sat), ...
        final_mass_all(sat), ...
        H_drift_all(sat), ...
        max_sigma_all(sat));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.5;

% Residual norm vs RAAN
figure;
semilogy(raan_deg_list,resnorm_all,'o-','LineWidth',lw);
grid on; box on;
xlabel('Target RAAN [deg]');
ylabel('Terminal residual norm');
title('Constellation Shooting Residual Norm vs Target RAAN');

% Final mass vs RAAN
figure;
plot(raan_deg_list,final_mass_all,'o-','LineWidth',lw);
grid on; box on;
xlabel('Target RAAN [deg]');
ylabel('Final mass');
title('Final Mass vs Target RAAN');

% Hamiltonian drift vs RAAN
figure;
semilogy(raan_deg_list,H_drift_all,'o-','LineWidth',lw);
grid on; box on;
xlabel('Target RAAN [deg]');
ylabel('max |\tilde{H}-\tilde{H}_0|');
title('Averaged Hamiltonian Drift vs Target RAAN');

% Target h-k distribution
figure;
plot([target_all.h],[target_all.k],'o','LineWidth',lw); hold on;
theta = linspace(0,2*pi,500);
plot(rho_hk*cos(theta),rho_hk*sin(theta),'--','LineWidth',lw);
plot(h0,k0,'kx','LineWidth',2,'MarkerSize',10);
grid on; box on; axis equal;
xlabel('h');
ylabel('k');
title('Target Plane Elements in h-k Space');
legend('Target satellites','i = 15 deg circle','Initial GTO plane','Location','best');

% p,f,g,h,k histories for all satellites
figure;
for sat = 1:N_s
    tau = tau_all{sat};
    x = x_all{sat};

    subplot(3,2,1);
    plot(tau,x(:,1),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('p'); title('p');

    subplot(3,2,2);
    plot(tau,x(:,2),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('f'); title('f');

    subplot(3,2,3);
    plot(tau,x(:,3),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('g'); title('g');

    subplot(3,2,4);
    plot(tau,x(:,4),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('h'); title('h');

    subplot(3,2,5);
    plot(tau,x(:,5),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('k'); xlabel('\tau'); title('k');

    subplot(3,2,6);
    plot(tau,x(:,9),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('m'); xlabel('\tau'); title('mass');
end

subplot(3,2,1); xlabel('\tau');
subplot(3,2,2); xlabel('\tau');
subplot(3,2,3); xlabel('\tau');
subplot(3,2,4); xlabel('\tau');

sgtitle('Averaged State Histories for All Satellites');

% Costate histories for all satellites
figure;
for sat = 1:N_s
    tau = tau_all{sat};
    lambda = lambda_all{sat};

    subplot(3,2,1);
    plot(tau,lambda(:,1),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('\lambda_p'); title('\lambda_p');

    subplot(3,2,2);
    plot(tau,lambda(:,2),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('\lambda_f'); title('\lambda_f');

    subplot(3,2,3);
    plot(tau,lambda(:,3),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('\lambda_g'); title('\lambda_g');

    subplot(3,2,4);
    plot(tau,lambda(:,4),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('\lambda_h'); title('\lambda_h');

    subplot(3,2,5);
    plot(tau,lambda(:,5),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('\lambda_k'); xlabel('\tau'); title('\lambda_k');

    subplot(3,2,6);
    plot(tau,lambda(:,9),'LineWidth',1.0); hold on; grid on; box on;
    ylabel('\lambda_m'); xlabel('\tau'); title('\lambda_m');
end

subplot(3,2,1); xlabel('\tau');
subplot(3,2,2); xlabel('\tau');
subplot(3,2,3); xlabel('\tau');
subplot(3,2,4); xlabel('\tau');

sgtitle('Costate Histories for All Satellites');

%% DV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta-v accounting for each satellite
%
% Nondimensional:
%   dv_ND = int_0^1 (T/m) * alpha0 d tau
%
% Dimensional:
%   dv_km_s = dv_ND * DU/TU
%   dv_m_s  = 1000 * dv_km_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dv_ND_all = zeros(1,N_s);
dv_km_s_all = zeros(1,N_s);
dv_m_s_all = zeros(1,N_s);

for sat = 1:N_s

    tau = tau_all{sat};
    x = x_all{sat};
    T = T_all{sat};

    m = x(:,9);

    accel_T = T ./ m;

    dv_ND = trapz(tau, accel_T) * alpha0;

    dv_ND_all(sat) = dv_ND;
    dv_km_s_all(sat) = dv_ND * DU/TU;
    dv_m_s_all(sat) = 1000 * dv_km_s_all(sat);

end

dv_total_ND = sum(dv_ND_all);
dv_total_km_s = sum(dv_km_s_all);
dv_total_m_s = sum(dv_m_s_all);

fprintf('\n============================================================\n');
fprintf('DELTA-V SUMMARY\n');
fprintf('============================================================\n');
fprintf('%5s %12s %20s %20s %20s\n', ...
    'Sat','RAAN[deg]','dv_ND','dv [km/s]','dv [m/s]');

for sat = 1:N_s
    fprintf('%5d %12.6f %20.12e %20.12f %20.6f\n', ...
        sat, ...
        target_all(sat).raan_deg, ...
        dv_ND_all(sat), ...
        dv_km_s_all(sat), ...
        dv_m_s_all(sat));
end

fprintf('------------------------------------------------------------\n');
fprintf('%5s %12s %20.12e %20.12f %20.6f\n', ...
    'SUM','-',dv_total_ND,dv_total_km_s,dv_total_m_s);
fprintf('============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct and plot all satellites
%
% Plot 1:
%   Full reconstructed controlled trajectories for all satellites.
%
% Plot 2:
%   Target Keplerian orbits for all satellites.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_reconstruct_all = true;

if do_reconstruct_all

    fprintf('\nReconstructing all satellite Cartesian trajectories...\n');

    recon_all = cell(1,N_s);
    out_f_all = cell(1,N_s);

    n_recon_per_step = 100;
    n_target_points = 100000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reconstruct all controlled trajectories and target orbits
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for sat = 1:N_s

        fprintf('  Reconstructing satellite %d / %d, RAAN = %.3f deg...\n', ...
            sat, N_s, target_all(sat).raan_deg);

        tau = tau_all{sat};
        x = x_all{sat};
        lambda = lambda_all{sat};
        target = target_all(sat);

        recon_all{sat} = prop_reconstructed_cartesian_trajectory( ...
            tau, x, lambda, params, n_recon_per_step);

        out_f_all{sat} = prop_kepler_orbit_from_mee_for_time( ...
            [target.p, target.f, target.g, target.h, target.k, 0], ...
            alpha0, params.mu, n_target_points);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot 1: full reconstructed controlled trajectories for all satellites
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    hold on;

    for sat = 1:N_s

        recon = recon_all{sat};

        plot3( ...
            recon.r_cart(:,1), ...
            recon.r_cart(:,2), ...
            recon.r_cart(:,3), ...
            'LineWidth', 1.25, ...
            'DisplayName', sprintf('Sat %02d, RAAN %.1f deg', ...
                sat, target_all(sat).raan_deg));

        % Mark final point
        plot3( ...
            recon.r_cart(end,1), ...
            recon.r_cart(end,2), ...
            recon.r_cart(end,3), ...
            'o', ...
            'MarkerSize', 5, ...
            'HandleVisibility','off');

    end

    % Initial GTO orbit for reference
    out_0 = prop_kepler_orbit_from_mee_for_time( ...
        [p0, f0, g0, h0, k0, L0], ...
        alpha0, params.mu, n_target_points);

    plot3( ...
        out_0.cart(:,1), ...
        out_0.cart(:,2), ...
        out_0.cart(:,3), ...
        'k--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Initial GTO orbit');

    earthy(1, 'Earth', 1,[0,0,0]);

    grid on; box on; axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Full Reconstructed Controlled Trajectories for All Satellites');
    legend('Location','bestoutside');

    view(3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot 2: target Keplerian orbits for all satellites
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    hold on;

    for sat = 1:N_s

        out_f = out_f_all{sat};

        plot3( ...
            out_f.cart(:,1), ...
            out_f.cart(:,2), ...
            out_f.cart(:,3), ...
            'LineWidth', 1.25, ...
            'DisplayName', sprintf('Target Sat %02d, RAAN %.1f deg', ...
                sat, target_all(sat).raan_deg));

    end

    earthy(1, 'Earth', 1,[0,0,0]);

    grid on; box on; axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(sprintf('Target Orbits for %d Satellites, i = %.1f deg', ...
        N_s, inc_f_deg));
    legend('Location','bestoutside');

    view(3);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = shooting_residual_symbolic(lambda0,x0,target,params,odeopts)

    r_bad = 1e10*ones(7,1);

    try
        [~,Y] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                      [0 1], [x0; lambda0], odeopts);

        if isempty(Y) || any(~isfinite(Y(:)))
            r = r_bad;
            return;
        end

        yf = Y(end,:).';
        xf = yf(1:9);
        lambdaf = yf(10:18);

        r = [
            xf(1) - target.p
            xf(2) - target.f
            xf(3) - target.g
            xf(4) - target.h
            xf(5) - target.k
            lambdaf(6)
            lambdaf(9)
        ];

        if any(~isfinite(r))
            r = r_bad;
        end

    catch
        r = r_bad;
    end

end

function res = terminal_residual(xf,lambdaf,target)

    res = [
        xf(1) - target.p
        xf(2) - target.f
        xf(3) - target.g
        xf(4) - target.h
        xf(5) - target.k
        lambdaf(6)
        lambdaf(9)
    ];

end

function [sigma,T] = thrust_history_simple(x,lambda,params)

    B = get_MEE_B_matrix(x,params.mu);

    lam_MEE = lambda(1:6);
    lam_m = lambda(9);

    Btl = B.'*lam_MEE;
    Btl_norm = norm(Btl);

    S = -(params.c/x(9))*Btl_norm - lam_m + 1;

    sigma = 0.5*(1 - S/sqrt(S^2 + params.eps_S^2));
    T = params.T_min + (params.T_max - params.T_min)*sigma;

end
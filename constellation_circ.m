%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% circularize_existing_constellation_from_workspace.m
%
% Assumes the eccentric RAAN-distributed constellation already exists in the
% workspace. This script does NOT clear variables.
%
% Required workspace variables:
%   x_all, lambda_all, target_all, params, DU, TU, N_s, odeopts, lsqopts
%
% Each satellite starts from its current final state:
%   x0_circ_sat = x_all{sat}(end,:).'
%
% Costate seed:
%   lambda0_seed = lambda_all{sat}(end,:).'
%
% Target:
%   circularize at current apogee radius:
%       p_circ = p_start/(1-e_start)
%       f = 0
%       g = 0
%       h,k unchanged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n============================================================\n');
fprintf('CIRCULARIZING EXISTING CONSTELLATION FROM WORKSPACE\n');
fprintf('============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_circ_steps = 15;          % continuation steps from eccentric to circular
use_previous_sat_seed = false;
% false: each satellite uses its own final costate from previous transfer
% true:  satellite j+1 uses optimized circularization costate from satellite j

alpha_circ_days = [];
% [] means keep the alpha already stored in each satellite final state.
% Otherwise set, for example:
% alpha_circ_days = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_circ_all = cell(1,N_s);
x_circ_all = cell(1,N_s);
lambda_circ_all = cell(1,N_s);
H_circ_all = cell(1,N_s);
sigma_circ_all = cell(1,N_s);
T_circ_all = cell(1,N_s);

lambda0_circ_all = zeros(9,N_s);
residual_circ_all = zeros(7,N_s);
resnorm_circ_all = zeros(1,N_s);

target_circ_all = struct([]);

final_mass_circ_all = zeros(1,N_s);
H_drift_circ_all = zeros(1,N_s);

dv_circ_ND_all = zeros(1,N_s);
dv_circ_km_s_all = zeros(1,N_s);
dv_circ_m_s_all = zeros(1,N_s);

lambda0_seed_global = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circularize each satellite
%
% Strategy:
%   Satellite 1:
%       use homotopy from eccentric target to circular target.
%
%   Satellites 2...N_s:
%       directly solve circular target using previous satellite's optimized
%       circularization costate as the initial guess.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sat = 1:N_s

    fprintf('\n============================================================\n');
    fprintf('CIRCULARIZATION SATELLITE %d / %d\n', sat, N_s);
    fprintf('============================================================\n');

    % Current final state and costate from previous constellation transfer
    x0_circ = x_all{sat}(end,:).';
    lambda0_from_previous_transfer = lambda_all{sat}(end,:).';

    % Optionally override circularization time of flight
    if ~isempty(alpha_circ_days)
        x0_circ(8) = alpha_circ_days*86400/TU;
    end

    % Current orbital shape
    p_start = x0_circ(1);
    f_start = x0_circ(2);
    g_start = x0_circ(3);
    h_start = x0_circ(4);
    k_start = x0_circ(5);

    e_start = sqrt(f_start^2 + g_start^2);

    if e_start >= 1
        error('Satellite %d has e >= 1. Cannot circularize elliptically.', sat);
    end

    % Circularize at current apogee radius
    p_circ = p_start/(1 - e_start);

    target_circ.p = p_circ;
    target_circ.f = 0;
    target_circ.g = 0;
    target_circ.h = h_start;
    target_circ.k = k_start;

    fprintf('Start p = %.15f\n', p_start);
    fprintf('Start f = %.15f\n', f_start);
    fprintf('Start g = %.15f\n', g_start);
    fprintf('Start e = %.15f\n', e_start);
    fprintf('Target circular p = %.15f\n', p_circ);
    fprintf('Target h = %.15f\n', target_circ.h);
    fprintf('Target k = %.15f\n', target_circ.k);
    fprintf('Circularization alpha = %.15f TU\n', x0_circ(8));
    fprintf('Circularization TOF = %.15f days\n', x0_circ(8)*TU/86400);

    if sat == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Satellite 1:
        % Homotopy from eccentric to circular.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fprintf('\nUsing circularization homotopy for satellite 1.\n');

        lambda0_work = lambda0_from_previous_transfer;

        beta_grid = linspace(0,1,N_circ_steps);

        for ib = 1:length(beta_grid)

            beta = beta_grid(ib);

            target.p = (1-beta)*p_start + beta*p_circ;
            target.f = (1-beta)*f_start;
            target.g = (1-beta)*g_start;
            target.h = h_start;
            target.k = k_start;

            fprintf('\n  Satellite %d circularization beta = %.6f\n', sat, beta);
            fprintf('    target.p = %.15f\n', target.p);
            fprintf('    target.f = %.15f\n', target.f);
            fprintf('    target.g = %.15f\n', target.g);
            fprintf('    target.h = %.15f\n', target.h);
            fprintf('    target.k = %.15f\n', target.k);

            lambda0_work = lsqnonlin( ...
                @(lambda0) shooting_residual_symbolic(lambda0,x0_circ,target,params,odeopts), ...
                lambda0_work, [], [], lsqopts);

            r_beta = shooting_residual_symbolic(lambda0_work,x0_circ,target,params,odeopts);

            fprintf('    residual norm = %.15e\n', norm(r_beta));

        end

    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Satellites 2...N_s:
        % Direct circularization using previous satellite's circularization
        % costate as the guess.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fprintf('\nUsing previous satellite circularization costate as seed.\n');

        lambda0_work = lambda0_seed_global;

        lambda0_work = lsqnonlin( ...
            @(lambda0) shooting_residual_symbolic(lambda0,x0_circ,target_circ,params,odeopts), ...
            lambda0_work, [], [], lsqopts);

        r_direct = shooting_residual_symbolic(lambda0_work,x0_circ,target_circ,params,odeopts);

        fprintf('Direct circularization residual norm = %.15e\n', norm(r_direct));

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final propagation for this satellite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [tau_circ,Y_circ] = ode45( ...
        @(tau,y) averaged_dynamics_symbolic(y,params), ...
        [0 1], [x0_circ; lambda0_work], odeopts);

    x_circ = Y_circ(:,1:9);
    lambda_circ = Y_circ(:,10:18);

    H_circ = zeros(length(tau_circ),1);
    sigma_circ = zeros(length(tau_circ),1);
    T_circ = zeros(length(tau_circ),1);

    for ii = 1:length(tau_circ)

        H_circ(ii) = calculate_averaged_hamiltonian_symbolic( ...
            x_circ(ii,:).',lambda_circ(ii,:).',params);

        [sigma_circ(ii),T_circ(ii)] = thrust_history_simple( ...
            x_circ(ii,:).',lambda_circ(ii,:).',params);

    end

    res_circ = terminal_residual( ...
        x_circ(end,:).',lambda_circ(end,:).',target_circ);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delta-v for circularization leg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m_hist = x_circ(:,9);
    accel_T = T_circ ./ m_hist;

    dv_ND = trapz(tau_circ, accel_T) * x0_circ(8);
    dv_km_s = dv_ND * DU/TU;
    dv_m_s = 1000*dv_km_s;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tau_circ_all{sat} = tau_circ;
    x_circ_all{sat} = x_circ;
    lambda_circ_all{sat} = lambda_circ;
    H_circ_all{sat} = H_circ;
    sigma_circ_all{sat} = sigma_circ;
    T_circ_all{sat} = T_circ;

    lambda0_circ_all(:,sat) = lambda0_work;
    residual_circ_all(:,sat) = res_circ;
    resnorm_circ_all(sat) = norm(res_circ);

    target_circ_all(sat).p = target_circ.p;
    target_circ_all(sat).f = target_circ.f;
    target_circ_all(sat).g = target_circ.g;
    target_circ_all(sat).h = target_circ.h;
    target_circ_all(sat).k = target_circ.k;

    if isfield(target_all,'raan_deg')
        target_circ_all(sat).raan_deg = target_all(sat).raan_deg;
    else
        target_circ_all(sat).raan_deg = NaN;
    end

    target_circ_all(sat).p_start = p_start;
    target_circ_all(sat).f_start = f_start;
    target_circ_all(sat).g_start = g_start;
    target_circ_all(sat).e_start = e_start;
    target_circ_all(sat).p_circ = p_circ;

    final_mass_circ_all(sat) = x_circ(end,9);
    H_drift_circ_all(sat) = max(abs(H_circ - H_circ(1)));

    dv_circ_ND_all(sat) = dv_ND;
    dv_circ_km_s_all(sat) = dv_km_s;
    dv_circ_m_s_all(sat) = dv_m_s;

    % This is the continuation seed for the next satellite
    lambda0_seed_global = lambda0_work;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print satellite summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n------------------------------------------------------------\n');
    fprintf('SATELLITE %d CIRCULARIZATION RESULT\n', sat);
    fprintf('------------------------------------------------------------\n');
    fprintf('terminal residual norm = %.15e\n', norm(res_circ));
    fprintf('Final p error = %.15e\n', x_circ(end,1)-target_circ.p);
    fprintf('Final f error = %.15e\n', x_circ(end,2)-target_circ.f);
    fprintf('Final g error = %.15e\n', x_circ(end,3)-target_circ.g);
    fprintf('Final h error = %.15e\n', x_circ(end,4)-target_circ.h);
    fprintf('Final k error = %.15e\n', x_circ(end,5)-target_circ.k);
    fprintf('lambda_L(tf) = %.15e\n', lambda_circ(end,6));
    fprintf('lambda_m(tf) = %.15e\n', lambda_circ(end,9));
    fprintf('Final mass = %.15f\n', x_circ(end,9));
    fprintf('max |H-H0| = %.15e\n', H_drift_circ_all(sat));
    fprintf('min sigma = %.15f\n', min(sigma_circ));
    fprintf('max sigma = %.15f\n', max(sigma_circ));
    fprintf('Delta-v = %.15e ND\n', dv_ND);
    fprintf('Delta-v = %.12f km/s\n', dv_km_s);
    fprintf('Delta-v = %.6f m/s\n', dv_m_s);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total circularization delta-v summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dv_circ_total_ND = sum(dv_circ_ND_all);
dv_circ_total_km_s = sum(dv_circ_km_s_all);
dv_circ_total_m_s = sum(dv_circ_m_s_all);

fprintf('\n============================================================\n');
fprintf('CIRCULARIZATION CONSTELLATION SUMMARY\n');
fprintf('============================================================\n');
fprintf('%5s %12s %20s %20s %20s %20s\n', ...
    'Sat','RAAN[deg]','resnorm','final mass','dv [m/s]','H drift');

for sat = 1:N_s

    if isfield(target_circ_all,'raan_deg')
        raan_print = target_circ_all(sat).raan_deg;
    else
        raan_print = NaN;
    end

    fprintf('%5d %12.6f %20.8e %20.12f %20.6f %20.8e\n', ...
        sat, ...
        raan_print, ...
        resnorm_circ_all(sat), ...
        final_mass_circ_all(sat), ...
        dv_circ_m_s_all(sat), ...
        H_drift_circ_all(sat));

end

fprintf('------------------------------------------------------------\n');
fprintf('%5s %12s %20s %20s %20.6f %20s\n', ...
    'SUM','-','-','-',dv_circ_total_m_s,'-');
fprintf('Total circularization delta-v = %.12f km/s\n', dv_circ_total_km_s);
fprintf('============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.5;

figure;
bar(dv_circ_m_s_all);
grid on; box on;
xlabel('Satellite index');
ylabel('\Delta v [m/s]');
title('Circularization \Delta v per Satellite');

figure;
semilogy(resnorm_circ_all,'o-','LineWidth',lw);
grid on; box on;
xlabel('Satellite index');
ylabel('Terminal residual norm');
title('Circularization Residual Norms');

figure;
plot(final_mass_circ_all,'o-','LineWidth',lw);
grid on; box on;
xlabel('Satellite index');
ylabel('Final mass');
title('Final Mass After Circularization');

figure;
hold on;
for sat = 1:N_s
    x_circ = x_circ_all{sat};
    tau_circ = tau_circ_all{sat};

    plot(tau_circ, sqrt(x_circ(:,2).^2 + x_circ(:,3).^2), ...
        'LineWidth',1.0);
end
grid on; box on;
xlabel('\tau');
ylabel('e = sqrt(f^2 + g^2)');
title('Eccentricity During Circularization');

%% Trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ALL SATELLITES:
%   1) full trajectory (transfer + circularization)
%   2) final orbits in inertial frame
%   3) final orbits in Earth-rotating frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_recon_per_step = 100;
n_orbit_pts = 4000;

omega_E_dim = 7.2921159e-5;   % rad/s
omega_E = omega_E_dim * TU;   % rad/TU

recon_transfer_all = cell(1,N_s);
recon_circ_all = cell(1,N_s);
orbit_inertial_all = cell(1,N_s);
orbit_rot_all = cell(1,N_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct everything first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sat = 1:N_s

    fprintf('Preparing plots for satellite %d / %d...\n', sat, N_s);

    %--------------------------------------------------
    % Reconstruct transfer leg
    %--------------------------------------------------
    recon_transfer_all{sat} = prop_reconstructed_cartesian_trajectory( ...
        tau_all{sat}, x_all{sat}, lambda_all{sat}, params, n_recon_per_step);

    %--------------------------------------------------
    % Reconstruct circularization leg
    %--------------------------------------------------
    recon_circ_all{sat} = prop_reconstructed_cartesian_trajectory( ...
        tau_circ_all{sat}, x_circ_all{sat}, lambda_circ_all{sat}, params, n_recon_per_step);

    %--------------------------------------------------
    % Final orbit in inertial frame
    % Use the actual final MEE state, including final longitude
    %--------------------------------------------------
    s_mee_f = x_circ_all{sat}(end,1:6).';

    p_f = s_mee_f(1);
    f_f = s_mee_f(2);
    g_f = s_mee_f(3);

    a_f = p_f/(1 - f_f^2 - g_f^2);
    T_orb = 2*pi*sqrt(a_f^3/params.mu);

    orbit_inertial_all{sat} = prop_kepler_orbit_from_mee_for_time( ...
        s_mee_f, T_orb, params.mu, n_orbit_pts);

    %--------------------------------------------------
    % Final orbit in rotating frame
    % Use the total elapsed mission time up to the end of circularization
    % as the Earth angle offset
    %--------------------------------------------------
    alpha_transfer = x_all{sat}(end,8);
    alpha_circ = x_circ_all{sat}(end,8);

    t_start_abs = alpha_transfer + alpha_circ;

    orbit_rot_all{sat} = rotate_cartesian_trajectory_to_earth_fixed( ...
        orbit_inertial_all{sat}.t, ...
        orbit_inertial_all{sat}.cart(:,1:3), ...
        omega_E, ...
        omega_E*t_start_abs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1: FULL TRAJECTORY OF ALL SATELLITES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

for sat = 1:N_s

    r_transfer = recon_transfer_all{sat}.r_cart;
    r_circ = recon_circ_all{sat}.r_cart;

    % Concatenate both legs, removing duplicate point
    r_full = [r_transfer; r_circ(2:end,:)];

    plot3( ...
        r_full(:,1), ...
        r_full(:,2), ...
        r_full(:,3), ...
        'LineWidth', 1.25, ...
        'DisplayName', sprintf('Sat %02d', sat));

end

earthy(1, 'Earth', 1,[0,0,0]);

grid on; box on; axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Full Trajectory of All Satellites (Transfer + Circularization)');
legend('Location','bestoutside');
view(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2: FINAL ORBITS IN INERTIAL FRAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

for sat = 1:N_s

    cart = orbit_inertial_all{sat}.cart;

    plot3( ...
        cart(:,1), ...
        cart(:,2), ...
        cart(:,3), ...
        'LineWidth', 1.25, ...
        'DisplayName', sprintf('Sat %02d', sat));

end

earthy(1, 'Earth', 1,[0,0,0]);

grid on; box on; axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Final Orbits in Inertial Frame');
legend('Location','bestoutside');
view(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3: FINAL ORBITS IN EARTH-ROTATING FRAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

for sat = 1:N_s

    r_rot = orbit_rot_all{sat}.r_rot;

    plot3( ...
        r_rot(:,1), ...
        r_rot(:,2), ...
        r_rot(:,3), ...
        'LineWidth', 1.25, ...
        'DisplayName', sprintf('Sat %02d', sat));

end

earthy(1, 'Earth', 1,[0,0,0]);

grid on; box on; axis equal;
xlabel('x_R');
ylabel('y_R');
zlabel('z_R');
title('Final Orbits in Earth-Rotating Frame');
legend('Location','bestoutside');
view(3);

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
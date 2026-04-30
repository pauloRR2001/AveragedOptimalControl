%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_paper_gto_geo_shooting_symbolic_lsqnonlin.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

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

% paper
% lambda0_guess = [
%    -2.321793758669676
%    -9.199953835732163
%     1.404963393724016
%     9.188784778243647
%    -1.546530651589456
%     0
%     0.000000000006312
%     0
%     0.074835258123558
% ];

% 2 hours run
% lambda0_guess = [
%   -1.828786152945906
%  -10.152640985483089
%    4.811173218653281
%   11.724568324458735
%   -0.898732369329855
%   -0.001893508407479
%    0.000006890251381
%                    0
%    0.038490613833436
% ];

% fast run with seps=1e-1, q=40 from fminunc
lambda0_guess = [
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

target.p = p0;
target.f = f0;
target.g = 0;
target.h = 0;
target.k = 0;

odeopts = odeset('RelTol',1e-12,'AbsTol',1e-12);

fprintf('\n========================================\n');
fprintf('PAPER GTO-TO-GTO-FLAT SHOOTING TEST (LSQNONLIN)\n');
fprintf('========================================\n');
fprintf('p0 = %.15f\n', x0(1));
fprintf('f0 = %.15f\n', x0(2));
fprintf('g0 = %.15f\n', x0(3));
fprintf('h0 = %.15f\n', x0(4));
fprintf('k0 = %.15f\n', x0(5));
fprintf('alpha0 = %.15f TU\n', x0(8));
fprintf('time of flight days = %.15f\n', x0(8)*TU/86400);
fprintf('m0 = %.15f kg\n', x0(9));
fprintf('params.c = %.15e\n', params.c);
fprintf('params.T_max = %.15e\n', params.T_max);
fprintf('params.T_max/m0 = %.15e\n', params.T_max/m0);
fprintf('params.J2 = %.15e\n', params.J2);
fprintf('params.q = %d\n', params.q);
fprintf('params.eps_S = %.15e\n', params.eps_S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixed propagation with seed costate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau_seed,Y_seed] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                          [0 1], [x0; lambda0_guess], odeopts);

x_seed = Y_seed(:,1:9);
lambda_seed = Y_seed(:,10:18);

H_seed = zeros(length(tau_seed),1);
sigma_seed = zeros(length(tau_seed),1);
T_seed = zeros(length(tau_seed),1);

for i = 1:length(tau_seed)
    H_seed(i) = calculate_averaged_hamiltonian_symbolic(x_seed(i,:).',lambda_seed(i,:).',params);
    [sigma_seed(i),T_seed(i)] = thrust_history_simple(x_seed(i,:).',lambda_seed(i,:).',params);
end

res_seed = shooting_residual_symbolic(lambda0_guess,x0,target,params,odeopts);

fprintf('\n========================================\n');
fprintf('FIXED PROPAGATION WITH SEED COSTATE\n');
fprintf('========================================\n');
fprintf('resnorm(seed) = %.15e\n', res_seed.'*res_seed);
fprintf('terminal residual norm = %.15e\n', norm(res_seed));
fprintf('Final p error = %.15e\n', x_seed(end,1)-target.p);
fprintf('Final f error = %.15e\n', x_seed(end,2)-target.f);
fprintf('Final g error = %.15e\n', x_seed(end,3)-target.g);
fprintf('Final h error = %.15e\n', x_seed(end,4)-target.h);
fprintf('Final k error = %.15e\n', x_seed(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda_seed(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda_seed(end,9));
fprintf('Final mass = %.15f\n', x_seed(end,9));
fprintf('max |H-H0| = %.15e\n', max(abs(H_seed-H_seed(1))));
fprintf('min sigma = %.15f\n', min(sigma_seed));
fprintf('max sigma = %.15f\n', max(sigma_seed));
fprintf('min T = %.15e\n', min(T_seed));
fprintf('max T = %.15e\n', max(T_seed));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lsqnonlin shooting solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lsqopts = optimoptions(@lsqnonlin, ...
    'Display','iter', ...
    'Algorithm','levenberg-marquardt', ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',5000, ...
    'StepTolerance',1e-12, ...
    'FunctionTolerance',1e-12, ...
    'OptimalityTolerance',1e-10);

lambda0_opt = lsqnonlin(@(lambda0) shooting_residual_symbolic(lambda0,x0,target,params,odeopts), ...
                        lambda0_guess, [], [], lsqopts);

[tau,Y] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                [0 1], [x0; lambda0_opt], odeopts);

x = Y(:,1:9);
lambda = Y(:,10:18);

H = zeros(length(tau),1);
sigma = zeros(length(tau),1);
T = zeros(length(tau),1);

for i = 1:length(tau)
    H(i) = calculate_averaged_hamiltonian_symbolic(x(i,:).',lambda(i,:).',params);
    [sigma(i),T(i)] = thrust_history_simple(x(i,:).',lambda(i,:).',params);
end

res_opt = terminal_residual(x(end,:).',lambda(end,:).',target);

%% Results

fprintf('\n========================================\n');
fprintf('LSQNONLIN CORRECTED SOLUTION\n');
fprintf('========================================\n');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.5;

figure;
plot(tau,x(:,1),'LineWidth',lw); hold on;
plot(tau,x(:,2),'LineWidth',lw);
plot(tau,x(:,3),'LineWidth',lw);
plot(tau,x(:,4),'LineWidth',lw);
plot(tau,x(:,5),'LineWidth',lw);
grid on; box on;
xlabel('$\tau$','Interpreter','latex');
ylabel('MEE states','Interpreter','latex');
title('Corrected GTO-to-GTO-Flat MEE State History','Interpreter','latex');
legend('p','f','g','h','k','Location','best','Interpreter','latex');

figure;
plot(tau,lambda(:,1),'LineWidth',lw); hold on;
plot(tau,lambda(:,2),'LineWidth',lw);
plot(tau,lambda(:,3),'LineWidth',lw);
plot(tau,lambda(:,4),'LineWidth',lw);
plot(tau,lambda(:,5),'LineWidth',lw);
plot(tau,lambda(:,6),'LineWidth',lw);
plot(tau,lambda(:,9),'LineWidth',lw);
grid on; box on;
xlabel('$\tau$','Interpreter','latex');
ylabel('Costates','Interpreter','latex');
title('Corrected GTO-to-GTO-Flat Costate History','Interpreter','latex');
legend('$\lambda_p$','$\lambda_f$','$\lambda_g$','$\lambda_h$','$\lambda_k$','$\lambda_L$','$\lambda_m$','Location','best'...
    ,'Interpreter','latex');

figure;
plot(tau,H-H(1),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('$\tilde H-\tilde H(0)$','Interpreter','latex');
title('Averaged Hamiltonian Error','Interpreter','latex');

figure;
plot(tau,sigma,'LineWidth',lw);
grid on; box on;
xlabel('$\tau$','Interpreter','latex');
ylabel('$\sigma$','Interpreter','latex');
title('Thrust Modulation','Interpreter','latex');

figure;
plot(tau,T,'LineWidth',lw);
grid on; box on;
xlabel('$\tau$','Interpreter','latex');
ylabel('T','Interpreter','latex');
title('Thrust Magnitude','Interpreter','latex');

figure;
plot(tau,x(:,9),'LineWidth',lw);
grid on; box on;
xlabel('$\tau$','Interpreter','latex');
ylabel('m','Interpreter','latex');
title('Mass History','Interpreter','latex');

%% reconstruct

recon = prop_reconstructed_cartesian_trajectory(tau, x, lambda, params, 100);

out_f = prop_kepler_orbit_from_mee_for_time([p0, f0, 0, 0, 0, 0], alpha0, params.mu, 100000);
out_0 = prop_kepler_orbit_from_mee_for_time([p0, f0, g0, h0, k0, L0], alpha0, params.mu, 100000);

%% final plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS USING recon AND out
%   recon.r_cart   = reconstructed transfer Cartesian history [N x 3]
%   recon.tau      = reconstructed transfer tau grid
%   recon.x        = reconstructed transfer state history [N x 9]
%
%   out.cart       = target orbit Cartesian history [M x 6] or [M x 3]
%   out.t          = target orbit time grid
%   out.kep        = target orbit Keplerian history [M x 6]
%
%   target         = struct with fields p,f,g,h,k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1: full trajectory and target orbit
figure;
plot3(recon.r_cart(:,1), recon.r_cart(:,2), recon.r_cart(:,3), 'b', 'LineWidth', 1.5); hold on;
plot3(out_f.cart(:,1),   out_f.cart(:,2),   out_f.cart(:,3),   'r', 'LineWidth', 1.5);
plot3(out_0.cart(:,1),   out_0.cart(:,2),   out_0.cart(:,3),   'k', 'LineWidth', 1.5);
earthy(1, 'Earth', 1,[0,0,0])
grid on; box on; axis equal;
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');
title('Transfer Trajectory and Target Orbit','Interpreter','latex');
legend('Controlled transfer','Target orbit','Initial Orbit','Location','best','Interpreter','latex');

% Figure 2: relative trajectory in LVLH-like frame built from target orbit
% Interpolate target Cartesian history onto reconstructed time grid
t_recon = recon.tau(:);
t_out   = out_f.t(:);

r_tgt_i = zeros(numel(t_recon),3);
v_tgt_i = zeros(numel(t_recon),3);

r_tgt_i(:,1) = interp1(t_out, out_f.cart(:,1), t_recon, 'pchip', 'extrap');
r_tgt_i(:,2) = interp1(t_out, out_f.cart(:,2), t_recon, 'pchip', 'extrap');
r_tgt_i(:,3) = interp1(t_out, out_f.cart(:,3), t_recon, 'pchip', 'extrap');

v_tgt_i(:,1) = interp1(t_out, out_f.cart(:,4), t_recon, 'pchip', 'extrap');
v_tgt_i(:,2) = interp1(t_out, out_f.cart(:,5), t_recon, 'pchip', 'extrap');
v_tgt_i(:,3) = interp1(t_out, out_f.cart(:,6), t_recon, 'pchip', 'extrap');

rho_lvlh = zeros(numel(t_recon),3);

for k = 1:numel(t_recon)

    rT = r_tgt_i(k,:).';
    vT = v_tgt_i(k,:).';
    rC = recon.r_cart(k,:).';

    er = rT / norm(rT);
    eh = cross(rT, vT);  eh = eh / norm(eh);
    et = cross(eh, er);

    C_I_to_RTN = [er.'; et.'; eh.'];

    rho_I = rC - rT;
    rho_lvlh(k,:) = (C_I_to_RTN * rho_I).';

end

figure;
plot3(rho_lvlh(:,1), rho_lvlh(:,2), rho_lvlh(:,3), 'k', 'LineWidth', 1.5);
grid on; box on; axis equal;
xlabel('$\Delta R$','Interpreter','latex');
ylabel('$\Delta T$','Interpreter','latex');
zlabel('$\Delta N$','Interpreter','latex');
title('Relative Trajectory with Respect to Target Orbit (RTN)','Interpreter','latex');
% optionally mark start/end
hold on;
plot3(rho_lvlh(1,1),   rho_lvlh(1,2),   rho_lvlh(1,3),   'go', 'MarkerFaceColor','g');
plot3(rho_lvlh(end,1), rho_lvlh(end,2), rho_lvlh(end,3), 'mo', 'MarkerFaceColor','m');
legend('Relative path','Start','End','Location','best','Interpreter','latex');

% Figure 3: state elements vs full transfer time with target final values

t_days = recon.tau * alpha0 * TU / 86400;

figure;

state_labels = {'p','f','g','h','k'};
target_vals = [target.p, target.f, target.g, target.h, target.k];

for j = 1:5
    subplot(3,2,j);
    plot(t_days, recon.x(:,j), 'b', 'LineWidth', 1.5); hold on;
    yline(target_vals(j), 'r--', 'LineWidth', 1.5);
    grid on; box on;
    xlabel('time [days]','Interpreter','latex');
    ylabel(state_labels{j},'Interpreter','latex');
    title([state_labels{j} ' vs time'],'Interpreter','latex');
    legend(state_labels{j}, 'target', 'Location', 'best','Interpreter','latex');
end

subplot(3,2,6);
plot(t_days, recon.x(:,9), 'b', 'LineWidth', 1.5);
grid on; box on;
xlabel('time [days]');
ylabel('m');
title('Mass vs time');
legend('m', 'Location', 'best');

%% modulation
lw = 1.5;

% figure;
% plot(recon.tau, recon.sigma, 'LineWidth', lw);
% grid on; box on;
% xlabel('\tau');
% ylabel('\sigma');
% title('Reconstructed Throttle History');
% 
% figure;
% plot(recon.tau, recon.T, 'LineWidth', lw);
% grid on; box on;
% xlabel('\tau');
% ylabel('T');
% title('Reconstructed Thrust Magnitude');

% figure;
% plot(recon.tau, recon.x(:,9), 'LineWidth', lw);
% grid on; box on;
% xlabel('\tau');
% ylabel('m');
% title('Reconstructed Mass History');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
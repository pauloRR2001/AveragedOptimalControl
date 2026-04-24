%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_lambda_consistency_along_trajectory.m
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

params.q = 20;
params.eps_S = 1e-1;
params.RE = 6378.1/DU;
params.J2 = 1.0826e-3;

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

lambda0 = [
   -1.3464
   -7.7480
    4.2458
   12.2186
    0.7624
   -0.0152
    0.0314
         0
    0.0393
];

odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 1: J2 = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.J2 = 0;

[tau0,Y0] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                  [0 1],[x0; lambda0],odeopts);

fprintf('\n========================================\n');
fprintf('CONSISTENCY ALONG TRAJECTORY: J2 = 0\n');
fprintf('========================================\n');

analyze_case(tau0,Y0,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 2: J2 = nominal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.J2 = 1.0826e-3;

[tau1,Y1] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                  [0 1],[x0; lambda0],odeopts);

fprintf('\n========================================\n');
fprintf('CONSISTENCY ALONG TRAJECTORY: J2 = nominal\n');
fprintf('========================================\n');

analyze_case(tau1,Y1,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function analyze_case(tau,Y,params)

    nsamp = 15;
    idx = round(linspace(1,length(tau),nsamp));

    err_mat = zeros(nsamp,9);
    Hvals = zeros(nsamp,1);

    fprintf('   tau          ||dlam||        e1             e2             e3             e4             e5             e6             e7             e8             e9\n');

    for k = 1:nsamp

        i = idx(k);

        y = Y(i,:).';
        x = y(1:9);
        lambda = y(10:18);

        ydot_prop = averaged_dynamics_symbolic(y,params);
        lambdadot_prop = ydot_prop(10:18);

        lambdadot_ham = lambdadot_from_sH(x,lambda,params);

        err = lambdadot_prop - lambdadot_ham;
        err_mat(k,:) = err.';
        Hvals(k) = calculate_averaged_hamiltonian_symbolic(x,lambda,params);

        fprintf('%8.5f   %12.5e  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n', ...
            tau(i), norm(err), err(1), err(2), err(3), err(4), err(5), err(6), err(7), err(8), err(9));

    end

    max_comp_err = max(abs(err_mat),[],1);

    fprintf('\nmax abs component mismatch along trajectory:\n');
    fprintf('e1 = %.15e\n', max_comp_err(1));
    fprintf('e2 = %.15e\n', max_comp_err(2));
    fprintf('e3 = %.15e\n', max_comp_err(3));
    fprintf('e4 = %.15e\n', max_comp_err(4));
    fprintf('e5 = %.15e\n', max_comp_err(5));
    fprintf('e6 = %.15e\n', max_comp_err(6));
    fprintf('e7 = %.15e\n', max_comp_err(7));
    fprintf('e8 = %.15e\n', max_comp_err(8));
    fprintf('e9 = %.15e\n', max_comp_err(9));

    fprintf('\nHamiltonian drift over sampled points:\n');
    fprintf('max |H-H0| = %.15e\n', max(abs(Hvals - Hvals(1))));

    figure;
    semilogy(tau(idx),sqrt(sum(err_mat.^2,2)),'o-','LineWidth',1.5);
    grid on; box on;
    xlabel('\tau');
    ylabel('||\Delta \lambda dot||');
    title('Costate Consistency Error Along Trajectory');

    figure;
    plot(tau(idx),err_mat(:,6),'o-','LineWidth',1.5); hold on;
    plot(tau(idx),err_mat(:,2),'o-','LineWidth',1.5);
    plot(tau(idx),err_mat(:,3),'o-','LineWidth',1.5);
    grid on; box on;
    xlabel('\tau');
    ylabel('component error');
    title('Key Costate Component Errors');
    legend('e_6 = \Delta \lambda_L','e_2 = \Delta \lambda_f','e_3 = \Delta \lambda_g','Location','best');

    figure;
    plot(tau(idx),Hvals - Hvals(1),'o-','LineWidth',1.5);
    grid on; box on;
    xlabel('\tau');
    ylabel('H-H(0)');
    title('Sampled Averaged Hamiltonian Drift');

end

function lambdadot_avg = lambdadot_from_sH(x,lambda,params)

    [z,w] = lgwt(params.q,-1,1);

    dHavg_dx = zeros(9,1);

    for i = 1:length(z)

        L = pi*z(i);

        x_i = x;
        x_i(6) = L;

        dsHdx_i = dsHdx_integrand_sym_fun( ...
            x_i, lambda, ...
            params.mu, params.c, params.T_min, params.T_max, ...
            params.eps_S, params.RE, params.J2);

        dHavg_dx = dHavg_dx + 0.5*w(i)*dsHdx_i;

    end

    lambdadot_avg = -dHavg_dx;

end
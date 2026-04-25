%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_avg_canonical_consistency.m
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
params.eps_S = 1e-2;

params.RE = 6378.1/DU;
params.J2 = 0;

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

y0 = [x0; lambda0];

odeopts = odeset('RelTol',1e-9,'AbsTol',1e-11);

fprintf('\n========================================\n');
fprintf('AVERAGED CANONICAL CONSISTENCY TEST\n');
fprintf('========================================\n');
fprintf('J2 = %.15e\n', params.J2);
fprintf('q = %d\n', params.q);
fprintf('eps_S = %.15e\n', params.eps_S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagate averaged system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau,Y] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), [0 1], y0, odeopts);

x = Y(:,1:9);
lambda = Y(:,10:18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hamiltonian drift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = zeros(length(tau),1);
for i = 1:length(tau)
    H(i) = calculate_averaged_hamiltonian_symbolic(x(i,:).',lambda(i,:).',params);
end

fprintf('\n========================================\n');
fprintf('AVERAGED PROPAGATION\n');
fprintf('========================================\n');
fprintf('max |H-H0| = %.15e\n', max(abs(H - H(1))));
fprintf('Final p = %.15f\n', x(end,1));
fprintf('Final f = %.15f\n', x(end,2));
fprintf('Final g = %.15f\n', x(end,3));
fprintf('Final h = %.15f\n', x(end,4));
fprintf('Final k = %.15f\n', x(end,5));
fprintf('Final L = %.15f\n', x(end,6));
fprintf('Final mass = %.15f\n', x(end,9));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% canonical consistency at sampled points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsamp = 15;
idx = round(linspace(1,length(tau),nsamp));

fprintf('\n========================================\n');
fprintf('POINTWISE CANONICAL CONSISTENCY\n');
fprintf('========================================\n');
fprintf('   tau          ||dx||          ||dlam||\n');

err_x = zeros(nsamp,1);
err_lam = zeros(nsamp,1);

for k = 1:nsamp

    i = idx(k);

    y = Y(i,:).';
    x_i = y(1:9);
    lambda_i = y(10:18);

    ydot_prop = averaged_dynamics_symbolic(y,params);
    xdot_prop = ydot_prop(1:9);
    lambdadot_prop = ydot_prop(10:18);

    xdot_ham = xdot_from_sH(x_i,lambda_i,params);
    lambdadot_ham = lambdadot_from_sH(x_i,lambda_i,params);

    err_x(k) = norm(xdot_prop - xdot_ham);
    err_lam(k) = norm(lambdadot_prop - lambdadot_ham);

    fprintf('%8.5f   %12.5e  %12.5e\n', tau(i), err_x(k), err_lam(k));

end

fprintf('\nmax ||xdot_prop - xdot_ham||   = %.15e\n', max(err_x));
fprintf('max ||lamdot_prop - lamdot_ham|| = %.15e\n', max(err_lam));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.5;

figure;
plot(tau,H - H(1),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('H_{avg} - H_{avg}(0)');
title('Averaged Hamiltonian Drift');

figure;
semilogy(tau,abs(H - H(1)) + 1e-30,'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('|H_{avg} - H_{avg}(0)|');
title('Absolute Averaged Hamiltonian Drift');

figure;
semilogy(tau(idx),err_x + 1e-30,'o-','LineWidth',lw); hold on;
semilogy(tau(idx),err_lam + 1e-30,'o-','LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('consistency error');
title('Pointwise Canonical Consistency');
legend('state error','costate error','Location','best');

function xdot_avg = xdot_from_sH(x,lambda,params)

    [z,w] = lgwt(params.q,-1,1);

    xdot_avg = zeros(9,1);

    for i = 1:length(z)

        L = pi*z(i);

        x_i = x;
        x_i(6) = L;

        dsHdlam_i = dsHdlam_integrand_sym_fun( ...
            x_i, lambda, ...
            params.mu, params.c, params.T_min, params.T_max, ...
            params.eps_S, params.RE, params.J2);

        xdot_avg = xdot_avg + 0.5*w(i)*dsHdlam_i;

    end

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
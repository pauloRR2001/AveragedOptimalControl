%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_avg_hamiltonian_gradient_check.m
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

[tau,Y] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), [0 1], y0, odeopts);

nsamp = 15;
idx = round(linspace(1,length(tau),nsamp));

fprintf('\n========================================\n');
fprintf('AVERAGED HAMILTONIAN GRADIENT CHECK\n');
fprintf('========================================\n');
fprintf('J2 = %.15e\n', params.J2);
fprintf('q = %d\n', params.q);
fprintf('eps_S = %.15e\n', params.eps_S);
fprintf('\n');
fprintf('   tau         ||dx-pH/plam||      ||dlam+pH/px||        |dH/dtau|          |H-H0|\n');

err_x = zeros(nsamp,1);
err_lam = zeros(nsamp,1);
err_Hdot = zeros(nsamp,1);
Herr = zeros(nsamp,1);

for k = 1:nsamp

    i = idx(k);

    y = Y(i,:).';
    x = y(1:9);
    lambda = y(10:18);

    ydot = averaged_dynamics_symbolic(y,params);
    xdot = ydot(1:9);
    lambdadot = ydot(10:18);

    dHdx = dHdx_avg_from_integrand(x,lambda,params);
    dHdlam = dHdlam_avg_from_integrand(x,lambda,params);

    err_x(k) = norm(xdot - dHdlam);
    err_lam(k) = norm(lambdadot + dHdx);

    H = calculate_averaged_hamiltonian_symbolic(x,lambda,params);
    if k == 1
        H0 = H;
    end
    Herr(k) = abs(H - H0);

    Hdot_chain = dHdx.'*xdot + dHdlam.'*lambdadot;
    err_Hdot(k) = abs(Hdot_chain);

    fprintf('%8.5f    %14.6e    %14.6e    %14.6e    %14.6e\n', ...
        tau(i), err_x(k), err_lam(k), err_Hdot(k), Herr(k));

end

fprintf('\nmax ||xdot - dH/dlambda|| = %.15e\n', max(err_x));
fprintf('max ||lambdadot + dH/dx|| = %.15e\n', max(err_lam));
fprintf('max |dH/dtau|             = %.15e\n', max(err_Hdot));
fprintf('max |H-H0| sampled        = %.15e\n', max(Herr));

figure;
semilogy(tau(idx),err_x + 1e-30,'o-','LineWidth',1.5); hold on;
semilogy(tau(idx),err_lam + 1e-30,'o-','LineWidth',1.5);
semilogy(tau(idx),err_Hdot + 1e-30,'o-','LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('error');
title('Averaged Hamiltonian Gradient Check');
legend('||xdot-dH/dlambda||','||lambdadot+dH/dx||','|dH/dtau|','Location','best');

figure;
plot(tau(idx),Herr,'o-','LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('|H-H(0)|');
title('Sampled Averaged Hamiltonian Drift');

function dHdx_avg = dHdx_avg_from_integrand(x,lambda,params)

[z,w] = lgwt(params.q,-1,1);

dHdx_avg = zeros(9,1);

for i = 1:length(z)

    L = pi*z(i);

    x_i = x;
    x_i(6) = L;

    dsHdx_i = dsHdx_integrand_sym_fun( ...
        x_i, lambda, ...
        params.mu, params.c, params.T_min, params.T_max, ...
        params.eps_S, params.RE, params.J2);

    dHdx_avg = dHdx_avg + 0.5*w(i)*dsHdx_i;

end

end

function dHdlam_avg = dHdlam_avg_from_integrand(x,lambda,params)

[z,w] = lgwt(params.q,-1,1);

dHdlam_avg = zeros(9,1);

for i = 1:length(z)

    L = pi*z(i);

    x_i = x;
    x_i(6) = L;

    dsHdlam_i = dsHdlam_integrand_sym_fun( ...
        x_i, lambda, ...
        params.mu, params.c, params.T_min, params.T_max, ...
        params.eps_S, params.RE, params.J2);

    dHdlam_avg = dHdlam_avg + 0.5*w(i)*dsHdlam_i;

end

end
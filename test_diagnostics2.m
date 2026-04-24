%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_symbolic_costate_consistency.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

DU = 6378;
TU = 806.78557;

params_base.mu = 1;
params_base.g0 = 9.80665e-3*(TU^2/DU);
params_base.Isp = 3100/TU;
params_base.c = params_base.g0*params_base.Isp;

params_base.T_min = 0;
params_base.T_max = 0.1*(TU^2/(DU*1000));

params_base.q = 20;
params_base.eps_S = 1e-1;

params_base.RE = 6378.1/DU;
params_base.J2 = 1.0826e-3;

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
   0.177539811977954
  -0.008044219729262
   4.592634444038394
  10.046597288646886
  -3.042766429481612
  -0.002942465450990
   0.085647894836700
  -0.085392440562130
   0.000746712757548
];

fprintf('\n========================================\n');
fprintf('SYMBOLIC COSTATE CONSISTENCY TEST\n');
fprintf('========================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 1: J2 = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = params_base;
params.J2 = 0;

ydot = averaged_dynamics_symbolic([x0; lambda0], params);
xdot_H = ydot(1:9);
lambdadot_H = ydot(10:18);

xdot_direct = direct_average_xdot_symbolic(x0, lambda0, params);
lambdadot_fd = numerical_lambdadot_from_Havg(x0, lambda0, params);

fprintf('\n========================================\n');
fprintf('CASE 1: J2 = 0\n');
fprintf('========================================\n');
fprintf('norm(xdot_H - xdot_direct) = %.15e\n', norm(xdot_H - xdot_direct));
fprintf('relative xdot error        = %.15e\n', norm(xdot_H - xdot_direct)/max(1,norm(xdot_direct)));
fprintf('norm(lambdadot_H - lambdadot_fd) = %.15e\n', norm(lambdadot_H - lambdadot_fd));
fprintf('relative lambdadot error        = %.15e\n', norm(lambdadot_H - lambdadot_fd)/max(1,norm(lambdadot_fd)));

disp('lambdadot_H (J2 = 0):');
disp(lambdadot_H.');

disp('lambdadot_fd (J2 = 0):');
disp(lambdadot_fd.');

disp('difference lambdadot_H - lambdadot_fd (J2 = 0):');
disp((lambdadot_H - lambdadot_fd).');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 2: J2 = nominal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = params_base;
params.J2 = 1.0826e-3;

ydot = averaged_dynamics_symbolic([x0; lambda0], params);
xdot_H = ydot(1:9);
lambdadot_H = ydot(10:18);

xdot_direct = direct_average_xdot_symbolic(x0, lambda0, params);
lambdadot_fd = numerical_lambdadot_from_Havg(x0, lambda0, params);

fprintf('\n========================================\n');
fprintf('CASE 2: J2 = nominal\n');
fprintf('========================================\n');
fprintf('norm(xdot_H - xdot_direct) = %.15e\n', norm(xdot_H - xdot_direct));
fprintf('relative xdot error        = %.15e\n', norm(xdot_H - xdot_direct)/max(1,norm(xdot_direct)));
fprintf('norm(lambdadot_H - lambdadot_fd) = %.15e\n', norm(lambdadot_H - lambdadot_fd));
fprintf('relative lambdadot error        = %.15e\n', norm(lambdadot_H - lambdadot_fd)/max(1,norm(lambdadot_fd)));

disp('lambdadot_H (J2 = nominal):');
disp(lambdadot_H.');

disp('lambdadot_fd (J2 = nominal):');
disp(lambdadot_fd.');

disp('difference lambdadot_H - lambdadot_fd (J2 = nominal):');
disp((lambdadot_H - lambdadot_fd).');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xdot_avg = direct_average_xdot_symbolic(x,lambda,params)

    [z,w] = lgwt(params.q,-1,1);
    xdot_avg = zeros(9,1);

    for i = 1:length(z)

        L = pi*z(i);

        x_i = x;
        x_i(6) = L;

        xdot_i = xdot_unavg_sym_fun( ...
            x_i, lambda, ...
            params.mu, params.c, params.T_min, params.T_max, ...
            params.eps_S, params.RE, params.J2);

        p = x_i(1);
        ff = x_i(2);
        g = x_i(3);

        q = 1 + ff*cos(L) + g*sin(L);
        r = p/q;
        dt_dL = r^2/sqrt(params.mu*p);
        a = p/(1 - ff^2 - g^2);
        P = 2*pi*sqrt(a^3/params.mu);
        s = (2*pi/P)*dt_dL;

        xdot_avg = xdot_avg + 0.5*w(i)*s*xdot_i;

    end

end

function lambdadot_fd = numerical_lambdadot_from_Havg(x,lambda,params)

    n = length(x);
    dHdx = zeros(n,1);
    dx = 1e-7;

    for j = 1:n
        xp = x;
        xm = x;

        xp(j) = xp(j) + dx;
        xm(j) = xm(j) - dx;

        Hp = calculate_averaged_hamiltonian_symbolic(xp,lambda,params);
        Hm = calculate_averaged_hamiltonian_symbolic(xm,lambda,params);

        dHdx(j) = (Hp - Hm)/(2*dx);
    end

    lambdadot_fd = -dHdx;

end
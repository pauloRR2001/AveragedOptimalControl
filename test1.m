%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test1.m
% Minimal brute-force averaged propagation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

DU = 6378;
TU = 806.78557;

params.mu = 1;
params.g0 = 9.80665e-3*(TU^2/DU);
params.Isp = 3100/TU;
params.c = params.g0*params.Isp;

params.T_min = 0;
params.T_max = 0.1*(TU^2/(100*DU*1000));

params.q = 20;
params.epsilon_S = 1e-1;
params.use_J2 = false;
params.RE = 6378.1/DU;
params.J2 = 1.0826e-3;

a0 = 24505/DU;
e0 = 0.725;
inc0 = deg2rad(28.5);

p0 = a0*(1 - e0^2);
f0 = e0;
g0 = 0;
h0 = tan(inc0/2);
k0 = 0;
L0 = 0;
t0 = 0;
alpha0 = 20*86400/TU;
m0 = 1;

x0 = [p0; f0; g0; h0; k0; L0; t0; alpha0; m0];

lambda0 = [
   -2.321793758669676
   -9.199953835732163
    1.404963393724016
    9.188784778243647
   -1.546530651589456
    0
    0.000000000006312
    0
    0.074835258123558
];

y0 = [x0; lambda0];

fprintf('p0 = %.15f\n', x0(1));
fprintf('f0 = %.15f\n', x0(2));
fprintf('g0 = %.15f\n', x0(3));
fprintf('h0 = %.15f\n', x0(4));
fprintf('k0 = %.15f\n', x0(5));
fprintf('alpha0 = %.15f\n', x0(8));
fprintf('m0 = %.15f\n', x0(9));
fprintf('params.c = %.15e\n', params.c);
fprintf('params.T_max = %.15e\n', params.T_max);

odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

[tau,Y] = ode45(@(tau,y) averaged_dynamics(y,params), [0 1], y0, odeopts);

x = Y(:,1:9);
lambda = Y(:,10:18);

H = zeros(length(tau),1);

for i = 1:length(tau)
    H(i) = calculate_averaged_hamiltonian(x(i,:).',lambda(i,:).',params);
end

sigma = zeros(length(tau),1);
T = zeros(length(tau),1);

for i = 1:length(tau)
    [sigma(i),T(i)] = thrust_history_simple(x(i,:).',lambda(i,:).',params);
end

fprintf('Final p = %.15f\n', x(end,1));
fprintf('Final f = %.15f\n', x(end,2));
fprintf('Final g = %.15f\n', x(end,3));
fprintf('Final h = %.15f\n', x(end,4));
fprintf('Final k = %.15f\n', x(end,5));
fprintf('Final L = %.15f\n', x(end,6));
fprintf('Final mass = %.15f\n', x(end,9));
fprintf('max abs(H-H0) = %.15e\n', max(abs(H-H(1))));

figure;
plot(tau,x(:,1),'LineWidth',1.5); hold on;
plot(tau,x(:,2),'LineWidth',1.5);
plot(tau,x(:,3),'LineWidth',1.5);
plot(tau,x(:,4),'LineWidth',1.5);
plot(tau,x(:,5),'LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('MEE states');
title('Brute-Force Averaged MEE States');
legend('p','f','g','h','k','Location','best');

figure;
plot(tau,lambda(:,1),'LineWidth',1.5); hold on;
plot(tau,lambda(:,2),'LineWidth',1.5);
plot(tau,lambda(:,3),'LineWidth',1.5);
plot(tau,lambda(:,4),'LineWidth',1.5);
plot(tau,lambda(:,5),'LineWidth',1.5);
plot(tau,lambda(:,6),'LineWidth',1.5);
plot(tau,lambda(:,9),'LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('Costates');
title('Brute-Force Averaged Costates');
legend('\lambda_p','\lambda_f','\lambda_g','\lambda_h','\lambda_k','\lambda_L','\lambda_m','Location','best');

figure;
plot(tau,H-H(1),'LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('\tilde H-\tilde H(0)');
title('Brute-Force Averaged Hamiltonian Error');

figure;
plot(tau,sigma,'LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('\sigma');
title('Thrust Modulation');

figure;
plot(tau,T,'LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('T');
title('Thrust Magnitude');

figure;
plot(tau,x(:,9),'LineWidth',1.5);
grid on; box on;
xlabel('\tau');
ylabel('m');
title('Mass History');

function [sigma,T] = thrust_history_simple(x,lambda,params)

    x_i = x;
    L = x_i(6);

    m = x_i(9);

    B = get_MEE_B_matrix(x_i,params.mu);

    lam_MEE = lambda(1:6);
    lam_m = lambda(9);

    Btl = B.'*lam_MEE;
    Btl_norm = norm(Btl);

    S = -(params.c/m)*Btl_norm - lam_m + 1;

    sigma = 0.5*(1 - S/sqrt(S^2 + params.epsilon_S^2));

    T = params.T_min + (params.T_max - params.T_min)*sigma;

end
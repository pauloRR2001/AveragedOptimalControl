%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 568 Final Project
% Main: Paper-Style Averaged MEE Primer-Vector Transfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

DU = 6378;
TU = 806.78557;

params.mu = 1;
params.g0 = 9.80665e-3*(TU^2/DU);
params.Isp = 3100/TU;
params.c = params.g0*params.Isp;
params.RE = 6378.1;
params.J2 = 1.0826e-3;

params.T_max = 0.1*(TU^2/(100*DU*1000)) / 1000;
params.q = 200;
params.epsilon_S = 1e-1;
params.T_min = 0;
params.use_J2 = true;

i0 = deg2rad(28.5);
Omega0 = 0;
omega0 = 0;
nu0 = 0;

DU = 6378;
TU = 806.78557;

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
alpha0 = 20*86400/TU / 20;
m0 = 1;

x0 = [p0; f0; g0; h0; k0; L0; t0; alpha0; m0];

target.p = p0;
target.f = f0;
target.g = g0;
target.h = 0;
target.k = 0;

lambda0_guess = [
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

odeopts = odeset('RelTol',1e-10,'AbsTol',1e-10);
optopts = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'OptimalityTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'MaxIterations',5000, ...
    'MaxFunctionEvaluations',50000);

lambda0_opt = fminunc(@(lambda0) shooting_cost(lambda0,x0,target,params,odeopts), ...
                      lambda0_guess,optopts);

[tau,Y] = ode45(@(tau,y) focused_averaged_dynamics(y,params), ...
                [0 1],[x0; lambda0_opt],odeopts);

x = Y(:,1:9);
lambda = Y(:,10:18);

H = zeros(length(tau),1);
for i = 1:length(tau)
    H(i) = calculate_focused_averaged_hamiltonian(x(i,:).',lambda(i,:).',params);
end

sigma = zeros(length(tau),1);
T = zeros(length(tau),1);

for i = 1:length(tau)
    [sigma(i),T(i)] = thrust_history(x(i,:).',lambda(i,:).',params);
end

fprintf('Optimized initial costates:\n');
disp(lambda0_opt);

fprintf('Final p error = %.15e\n', x(end,1)-target.p);
fprintf('Final f error = %.15e\n', x(end,2)-target.f);
fprintf('Final g error = %.15e\n', x(end,3)-target.g);
fprintf('Final h error = %.15e\n', x(end,4)-target.h);
fprintf('Final k error = %.15e\n', x(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda(end,9));
fprintf('Final mass = %.15f\n', x(end,9));
fprintf('Hamiltonian RMS error = %.15e\n', rms((H-H(1))/H(1)));

lw = 1.8;
fs = 14;

figure;
plot(tau,x(:,1),'LineWidth',lw); hold on;
plot(tau,x(:,2),'LineWidth',lw);
plot(tau,x(:,3),'LineWidth',lw);
plot(tau,x(:,4),'LineWidth',lw);
plot(tau,x(:,5),'LineWidth',lw);
grid on; box on;
xlabel('\tau','FontSize',fs);
ylabel('MEE States','FontSize',fs);
title('Averaged MEE State History','FontSize',fs);
legend('p','f','g','h','k','Location','best');
set(gca,'FontSize',fs);

figure;
plot(tau,lambda(:,1),'LineWidth',lw); hold on;
plot(tau,lambda(:,2),'LineWidth',lw);
plot(tau,lambda(:,3),'LineWidth',lw);
plot(tau,lambda(:,4),'LineWidth',lw);
plot(tau,lambda(:,5),'LineWidth',lw);
plot(tau,lambda(:,6),'LineWidth',lw);
plot(tau,lambda(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau','FontSize',fs);
ylabel('Costates','FontSize',fs);
title('Averaged Costate History','FontSize',fs);
legend('\lambda_p','\lambda_f','\lambda_g','\lambda_h','\lambda_k','\lambda_L','\lambda_m','Location','best');
set(gca,'FontSize',fs);

figure;
plot(tau,T,'LineWidth',lw);
grid on; box on;
xlabel('\tau','FontSize',fs);
ylabel('Thrust','FontSize',fs);
title('Smoothed Thrust History','FontSize',fs);
set(gca,'FontSize',fs);

figure;
plot(tau,sigma,'LineWidth',lw);
grid on; box on;
xlabel('\tau','FontSize',fs);
ylabel('\sigma','FontSize',fs);
title('Thrust Modulation History','FontSize',fs);
set(gca,'FontSize',fs);

figure;
plot(tau,H-H(1),'LineWidth',lw);
grid on; box on;
xlabel('\tau','FontSize',fs);
ylabel('\tilde{H}-\tilde{H}(0)','FontSize',fs);
title('Averaged Hamiltonian Error','FontSize',fs);
set(gca,'FontSize',fs);

figure;
plot(tau,x(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau','FontSize',fs);
ylabel('m','FontSize',fs);
title('Mass History','FontSize',fs);
set(gca,'FontSize',fs);

function J = shooting_cost(lambda0,x0,target,params,odeopts)

[tau,Y] = ode45(@(tau,y) focused_averaged_dynamics(y,params), ...
                [0 1],[x0; lambda0],odeopts);

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

J = r.'*r;

end

function [sigma,T] = thrust_history(x,lambda,params)

m = x(9);

B = get_MEE_B_matrix(x,params.mu);
lam_MEE = lambda(1:6);

Btl = B.'*lam_MEE;
Btl_norm = norm(Btl);

S = -(params.c/m)*Btl_norm - lambda(9) + 1;

sigma = 0.5*(1 - S/sqrt(S^2 + params.epsilon_S^2));

T = params.T_min + (params.T_max - params.T_min)*sigma;

end
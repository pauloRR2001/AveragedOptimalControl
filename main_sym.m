%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_paper_gto_geo_shooting_symbolic.m
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
lambda0_guess =[
  -1.828786152945906
 -10.152640985483089
   4.811173218653281
  11.724568324458735
  -0.898732369329855
  -0.001893508407479
   0.000006890251381
                   0
   0.038490613833436
];

target.p = p0;
target.f = f0;
target.g = 0;
target.h = 0;
target.k = 0;

odeopts = odeset('RelTol',1e-12,'AbsTol',1e-12);

fprintf('\n========================================\n');
fprintf('PAPER GTO-TO-GTO-FLAT SHOOTING TEST\n');
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
% fixed propagation with appendix costate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau0,Y0] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                  [0 1], [x0; lambda0_guess], odeopts);

x_fixed = Y0(:,1:9);
lambda_fixed = Y0(:,10:18);

H_fixed = zeros(length(tau0),1);
sigma_fixed = zeros(length(tau0),1);
T_fixed = zeros(length(tau0),1);

for i = 1:length(tau0)
    H_fixed(i) = calculate_averaged_hamiltonian_symbolic(x_fixed(i,:).',lambda_fixed(i,:).',params);
    [sigma_fixed(i),T_fixed(i)] = thrust_history_simple(x_fixed(i,:).',lambda_fixed(i,:).',params);
end

J_fixed = shooting_cost(lambda0_guess,x0,target,params,odeopts);
res_fixed = terminal_residual(x_fixed(end,:).',lambda_fixed(end,:).',target);

fprintf('\n========================================\n');
fprintf('FIXED PROPAGATION WITH APPENDIX COSTATE\n');
fprintf('========================================\n');
fprintf('Z(lambda0_guess) = %.15e\n', J_fixed);
fprintf('terminal residual norm = %.15e\n', norm(res_fixed));
fprintf('Final p error = %.15e\n', x_fixed(end,1)-target.p);
fprintf('Final f error = %.15e\n', x_fixed(end,2)-target.f);
fprintf('Final g error = %.15e\n', x_fixed(end,3)-target.g);
fprintf('Final h error = %.15e\n', x_fixed(end,4)-target.h);
fprintf('Final k error = %.15e\n', x_fixed(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda_fixed(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda_fixed(end,9));
fprintf('Final mass = %.15f\n', x_fixed(end,9));
fprintf('max |H-H0| = %.15e\n', max(abs(H_fixed-H_fixed(1))));
fprintf('min sigma = %.15f\n', min(sigma_fixed));
fprintf('max sigma = %.15f\n', max(sigma_fixed));
fprintf('min T = %.15e\n', min(T_fixed));
fprintf('max T = %.15e\n', max(T_fixed));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shooting correction with paper objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fminopts = optimoptions(@fminunc, ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',5000, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-10, ...
    'FunctionTolerance',1e-12);

opts = optimoptions(@fminunc, ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',400, ...
    'MaxFunctionEvaluations',20000, ...
    'StepTolerance',1e-14, ...
    'OptimalityTolerance',1e-8, ...
    'FunctionTolerance',1e-14);

lambda0_opt = fminunc(@(lambda0) shooting_cost(lambda0,x0,target,params,odeopts), ...
                      lambda0_guess, fminopts);

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

J_opt = shooting_cost(lambda0_opt,x0,target,params,odeopts);
res_opt = terminal_residual(x(end,:).',lambda(end,:).',target);

fprintf('\n========================================\n');
fprintf('CORRECTED SHOOTING SOLUTION\n');
fprintf('========================================\n');
fprintf('Optimized initial costates:\n');
disp(lambda0_opt);
fprintf('Z(lambda0_opt) = %.15e\n', J_opt);
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
% mean Cartesian trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_cart = zeros(length(tau),3);

for i = 1:length(tau)
    s_mee = x(i,1:6).';
    s_kep = MEE2KEP(s_mee);
    s_cart = KEP2CART(s_kep,params.mu);
    r_cart(i,:) = s_cart(1:3).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.5;

figure;
plot(tau,x(:,1),'LineWidth',lw); hold on;
plot(tau,x(:,2),'LineWidth',lw);
plot(tau,x(:,3),'LineWidth',lw);
plot(tau,x(:,4),'LineWidth',lw);
plot(tau,x(:,5),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('MEE states');
title('Corrected GTO-to-GTO-Flat MEE State History');
legend('p','f','g','h','k','Location','best');

figure;
plot(tau,lambda(:,1),'LineWidth',lw); hold on;
plot(tau,lambda(:,2),'LineWidth',lw);
plot(tau,lambda(:,3),'LineWidth',lw);
plot(tau,lambda(:,4),'LineWidth',lw);
plot(tau,lambda(:,5),'LineWidth',lw);
plot(tau,lambda(:,6),'LineWidth',lw);
plot(tau,lambda(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('Costates');
title('Corrected GTO-to-GTO-Flat Costate History');
legend('\lambda_p','\lambda_f','\lambda_g','\lambda_h','\lambda_k','\lambda_L','\lambda_m','Location','best');

figure;
plot(tau,H-H(1),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('\tilde H-\tilde H(0)');
title('Averaged Hamiltonian Error');

figure;
plot(tau,sigma,'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('\sigma');
title('Thrust Modulation');

figure;
plot(tau,T,'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('T');
title('Thrust Magnitude');

figure;
plot(tau,x(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('m');
title('Mass History');

figure;
plot3(r_cart(:,1),r_cart(:,2),r_cart(:,3),'LineWidth',lw);
grid on; box on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Mean Cartesian Trajectory from Averaged MEE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = shooting_cost(lambda0,x0,target,params,odeopts)

    [~,Y] = ode45(@(tau,y) averaged_dynamics_symbolic(y,params), ...
                  [0 1], [x0; lambda0], odeopts);

    yf = Y(end,:).';
    xf = yf(1:9);
    lambdaf = yf(10:18);

    e = [
        xf(1) - target.p
        xf(2) - target.f
        xf(3) - target.g
        xf(4) - target.h
        xf(5) - target.k
    ];

    J = e.'*e + lambdaf(9)^2 + lambdaf(6)^2;

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
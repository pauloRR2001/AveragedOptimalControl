%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_dimensional_paper_scaling.m
% Diagnostics using paper-style mass scaling: m0 = 100 kg
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
m0 = 100;

x0 = [p0; f0; g0; h0; k0; L0; t0; alpha0; m0];

lambda_paper = [
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

% mid correction lambda
lambda_current = [
  -2.174262969144880
  -9.247345880446717
   1.437573416431035
   9.315768483082604
  -1.523339873735115
   0.038186148741821
  -0.066617183878231
                   0
   0.260025257864533
];

target.p = p0;
target.f = f0;
target.g = g0;
target.h = 0;
target.k = 0;

odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

fprintf('\n========================================\n');
fprintf('BASELINE VALUES\n');
fprintf('========================================\n');
fprintf('p0 = %.15f\n', x0(1));
fprintf('f0 = %.15f\n', x0(2));
fprintf('g0 = %.15f\n', x0(3));
fprintf('h0 = %.15f\n', x0(4));
fprintf('k0 = %.15f\n', x0(5));
fprintf('alpha0 = %.15f\n', x0(8));
fprintf('m0 = %.15f\n', x0(9));
fprintf('params.c = %.15e\n', params.c);
fprintf('params.T_max = %.15e\n', params.T_max);
fprintf('params.T_max/m0 = %.15e\n', params.T_max/m0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: zero-thrust two-body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_zero = params;
params_zero.T_min = 0;
params_zero.T_max = 0;
params_zero.use_J2 = false;

[tau1,Y1] = ode45(@(tau,y) averaged_dynamics(y,params_zero), ...
                  [0 1],[x0; lambda_paper],odeopts);

x1 = Y1(:,1:9);
lambda1 = Y1(:,10:18);

H1 = zeros(length(tau1),1);
for i = 1:length(tau1)
    H1(i) = calculate_averaged_hamiltonian(x1(i,:).',lambda1(i,:).',params_zero);
end

fprintf('\n========================================\n');
fprintf('TEST 1: ZERO-THRUST TWO-BODY SANITY\n');
fprintf('========================================\n');
fprintf('max |p-p0| = %.15e\n', max(abs(x1(:,1)-x0(1))));
fprintf('max |f-f0| = %.15e\n', max(abs(x1(:,2)-x0(2))));
fprintf('max |g-g0| = %.15e\n', max(abs(x1(:,3)-x0(3))));
fprintf('max |h-h0| = %.15e\n', max(abs(x1(:,4)-x0(4))));
fprintf('max |k-k0| = %.15e\n', max(abs(x1(:,5)-x0(5))));
fprintf('max |m-m0| = %.15e\n', max(abs(x1(:,9)-x0(9))));
fprintf('max |H-H0| = %.15e\n', max(abs(H1-H1(1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: Hamiltonian xdot vs direct average xdot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ydot_H = averaged_dynamics([x0; lambda_paper],params);
xdot_H = ydot_H(1:9);

xdot_direct = direct_average_xdot(x0,lambda_paper,params);

fprintf('\n========================================\n');
fprintf('TEST 2: HAMILTONIAN XDOT VS DIRECT AVG XDOT\n');
fprintf('========================================\n');
fprintf('norm(xdot_H - xdot_direct) = %.15e\n', norm(xdot_H - xdot_direct));
fprintf('relative error = %.15e\n', norm(xdot_H - xdot_direct)/max(1,norm(xdot_direct)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: fixed paper costate propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau3,Y3] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                  [0 1],[x0; lambda_paper],odeopts);

x3 = Y3(:,1:9);
lambda3 = Y3(:,10:18);

H3 = zeros(length(tau3),1);
sigma3 = zeros(length(tau3),1);
T3 = zeros(length(tau3),1);

for i = 1:length(tau3)
    H3(i) = calculate_averaged_hamiltonian(x3(i,:).',lambda3(i,:).',params);
    [sigma3(i),T3(i)] = thrust_history_simple(x3(i,:).',lambda3(i,:).',params);
end

res3 = terminal_residual(x3(end,:).',lambda3(end,:).',target);

fprintf('\n========================================\n');
fprintf('TEST 3: FIXED PAPER COSTATE PROPAGATION\n');
fprintf('========================================\n');
fprintf('terminal residual norm = %.15e\n', norm(res3));
fprintf('Final p error = %.15e\n', x3(end,1)-target.p);
fprintf('Final f error = %.15e\n', x3(end,2)-target.f);
fprintf('Final g error = %.15e\n', x3(end,3)-target.g);
fprintf('Final h error = %.15e\n', x3(end,4)-target.h);
fprintf('Final k error = %.15e\n', x3(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda3(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda3(end,9));
fprintf('Final mass = %.15f\n', x3(end,9));
fprintf('max |H-H0| = %.15e\n', max(abs(H3-H3(1))));
fprintf('min sigma = %.15f\n', min(sigma3));
fprintf('max sigma = %.15f\n', max(sigma3));

%% 3.5
lsqopts = optimoptions(@lsqnonlin, ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxIterations',100, ...
    'MaxFunctionEvaluations',2000);

lambda_opt = lsqnonlin(@(lambda0) shooting_residual(lambda0,x0,target,params,odeopts), ...
                       lambda_current, [], [], lsqopts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: optimize costates with paper scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optopts = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'OptimalityTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',5000);

lambda_opt = fminunc(@(lambda0) shooting_cost(lambda0,x0,target,params,odeopts), ...
                     lambda_paper,optopts);

[tau4,Y4] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                  [0 1],[x0; lambda_opt],odeopts);

x4 = Y4(:,1:9);
lambda4 = Y4(:,10:18);

H4 = zeros(length(tau4),1);
sigma4 = zeros(length(tau4),1);
T4 = zeros(length(tau4),1);

for i = 1:length(tau4)
    H4(i) = calculate_averaged_hamiltonian(x4(i,:).',lambda4(i,:).',params);
    [sigma4(i),T4(i)] = thrust_history_simple(x4(i,:).',lambda4(i,:).',params);
end

res4 = terminal_residual(x4(end,:).',lambda4(end,:).',target);

fprintf('\n========================================\n');
fprintf('TEST 4: OPTIMIZED COSTATES, PAPER MASS SCALING\n');
fprintf('========================================\n');
fprintf('Optimized costates:\n');
disp(lambda_opt);
fprintf('terminal residual norm = %.15e\n', norm(res4));
fprintf('Final p error = %.15e\n', x4(end,1)-target.p);
fprintf('Final f error = %.15e\n', x4(end,2)-target.f);
fprintf('Final g error = %.15e\n', x4(end,3)-target.g);
fprintf('Final h error = %.15e\n', x4(end,4)-target.h);
fprintf('Final k error = %.15e\n', x4(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda4(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda4(end,9));
fprintf('Final mass = %.15f\n', x4(end,9));
fprintf('max |H-H0| = %.15e\n', max(abs(H4-H4(1))));
fprintf('min sigma = %.15f\n', min(sigma4));
fprintf('max sigma = %.15f\n', max(sigma4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 5: continuation in epsilon_S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_list = [1e-1 5e-2 1e-2];

lambda_guess = lambda_opt;

fprintf('\n========================================\n');
fprintf('TEST 5: EPSILON CONTINUATION\n');
fprintf('========================================\n');

eps_res = zeros(length(eps_list),1);

for j = 1:length(eps_list)

    params_eps = params;
    params_eps.epsilon_S = eps_list(j);

    lambda_guess = fminunc(@(lambda0) shooting_cost(lambda0,x0,target,params_eps,odeopts), ...
                           lambda_guess,optopts);

    [tau_eps,Y_eps] = ode45(@(tau,y) averaged_dynamics(y,params_eps), ...
                            [0 1],[x0; lambda_guess],odeopts);

    x_eps = Y_eps(:,1:9);
    lambda_eps = Y_eps(:,10:18);

    res_eps = terminal_residual(x_eps(end,:).',lambda_eps(end,:).',target);

    eps_res(j) = norm(res_eps);

    fprintf('epsilon_S = %.3e | resnorm = %.15e | final mass = %.15f\n', ...
            eps_list(j), eps_res(j), x_eps(end,9));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.5;

figure;
plot(tau3,x3(:,1),'LineWidth',lw); hold on;
plot(tau3,x3(:,2),'LineWidth',lw);
plot(tau3,x3(:,3),'LineWidth',lw);
plot(tau3,x3(:,4),'LineWidth',lw);
plot(tau3,x3(:,5),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('MEE states');
title('Test 3: Fixed Paper Costate State History');
legend('p','f','g','h','k','Location','best');

figure;
plot(tau4,x4(:,1),'LineWidth',lw); hold on;
plot(tau4,x4(:,2),'LineWidth',lw);
plot(tau4,x4(:,3),'LineWidth',lw);
plot(tau4,x4(:,4),'LineWidth',lw);
plot(tau4,x4(:,5),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('MEE states');
title('Test 4: Optimized State History');
legend('p','f','g','h','k','Location','best');

figure;
plot(tau4,lambda4(:,1),'LineWidth',lw); hold on;
plot(tau4,lambda4(:,2),'LineWidth',lw);
plot(tau4,lambda4(:,3),'LineWidth',lw);
plot(tau4,lambda4(:,4),'LineWidth',lw);
plot(tau4,lambda4(:,5),'LineWidth',lw);
plot(tau4,lambda4(:,6),'LineWidth',lw);
plot(tau4,lambda4(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('Costates');
title('Test 4: Optimized Costate History');
legend('\lambda_p','\lambda_f','\lambda_g','\lambda_h','\lambda_k','\lambda_L','\lambda_m','Location','best');

figure;
plot(tau4,H4-H4(1),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('\tilde H-\tilde H(0)');
title('Test 4: Hamiltonian Error');

figure;
plot(tau4,sigma4,'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('\sigma');
title('Test 4: Thrust Modulation');

figure;
plot(tau4,x4(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('m');
title('Test 4: Mass History');

figure;
semilogy(eps_list,eps_res,'o-','LineWidth',lw);
grid on; box on;
xlabel('\epsilon_S');
ylabel('terminal residual norm');
title('Test 5: Epsilon Continuation Residuals');
set(gca,'XDir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = shooting_cost(lambda0,x0,target,params,odeopts)

    [~,Y] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                  [0 1],[x0; lambda0],odeopts);

    yf = Y(end,:).';
    xf = yf(1:9);
    lambdaf = yf(10:18);

    res = terminal_residual(xf,lambdaf,target);

    W = diag([1 1 1 10 10 1 1]);

    J = res.'*W*res;

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

    sigma = 0.5*(1 - S/sqrt(S^2 + params.epsilon_S^2));

    T = params.T_min + (params.T_max - params.T_min)*sigma;

end

function xdot_avg = direct_average_xdot(x,lambda,params)

    [z,w] = lgwt(params.q,-1,1);

    xdot_avg = zeros(9,1);

    for i = 1:length(z)

        L = pi*z(i);

        x_i = x;
        x_i(6) = L;

        alpha = x_i(8);
        m = x_i(9);

        A = get_MEE_A_matrix(x_i,params.mu);
        B = get_MEE_B_matrix(x_i,params.mu);

        if params.use_J2
            gamma_RTN = J2pert(x_i,params.mu,params.RE,params.J2);
        else
            gamma_RTN = [0;0;0];
        end

        lam_MEE = lambda(1:6);
        lam_m = lambda(9);

        Btl = B.'*lam_MEE;
        Btl_norm = norm(Btl);

        if Btl_norm < 1e-12
            u_hat = [0;0;0];
        else
            u_hat = -Btl/Btl_norm;
        end

        S = -(params.c/m)*Btl_norm - lam_m + 1;
        sigma = 0.5*(1 - S/sqrt(S^2 + params.epsilon_S^2));
        T = params.T_min + (params.T_max - params.T_min)*sigma;

        f_MEE = alpha*(A + B*(u_hat*(T/m) + gamma_RTN));
        f_t = alpha;
        f_alpha = 0;
        f_m = -alpha*T/params.c;

        f_full = [f_MEE; f_t; f_alpha; f_m];

        p = x_i(1);
        f = x_i(2);
        g = x_i(3);

        q = 1 + f*cos(L) + g*sin(L);
        r = p/q;
        dt_dL = r^2/sqrt(params.mu*p);
        a = p/(1 - f^2 - g^2);
        P = 2*pi*sqrt(a^3/params.mu);
        s = (2*pi/P)*dt_dL;

        xdot_avg = xdot_avg + 0.5*w(i)*s*f_full;

    end

end

function res = shooting_residual(lambda0,x0,target,params,odeopts)

    [~,Y] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                  [0 1],[x0; lambda0],odeopts);

    yf = Y(end,:).';
    xf = yf(1:9);
    lambdaf = yf(10:18);

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
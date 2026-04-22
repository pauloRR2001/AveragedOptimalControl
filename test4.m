%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_gto_paper.m
% GTO-to-GEO flat test using paper initial conditions and paper costates
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
params.use_J2 = true;
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

lambda0_paper_avg = [
   -2.321793758669676
   -9.199538357321630
    1.404633937240160
    9.188477824364700
   -1.546306515894560
    0
    0.000000000006312
    0
    0.074835258123558
];

target.p = p0;
target.f = f0;
target.g = g0;
target.h = 0;
target.k = 0;

odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

fprintf('\n========================================\n');
fprintf('PAPER GTO-TO-GEO AVG TEST\n');
fprintf('========================================\n');
fprintf('p0 = %.15f\n', x0(1));
fprintf('f0 = %.15f\n', x0(2));
fprintf('g0 = %.15f\n', x0(3));
fprintf('h0 = %.15f\n', x0(4));
fprintf('k0 = %.15f\n', x0(5));
fprintf('alpha0 = %.15f TU\n', x0(8));
fprintf('time of flight = %.15f days\n', x0(8)*TU/86400);
fprintf('m0 = %.15f kg\n', x0(9));
fprintf('params.c = %.15e\n', params.c);
fprintf('params.T_max = %.15e\n', params.T_max);
fprintf('params.T_max/m0 = %.15e\n', params.T_max/m0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1: zero-thrust sanity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_zero = params;
params_zero.T_min = 0;
params_zero.T_max = 0;
params_zero.use_J2 = false;

[tau_z,Y_z] = ode45(@(tau,y) averaged_dynamics(y,params_zero), ...
                    [0 1],[x0; lambda0_paper_avg],odeopts);

x_z = Y_z(:,1:9);
lambda_z = Y_z(:,10:18);

H_z = zeros(length(tau_z),1);
for i = 1:length(tau_z)
    H_z(i) = calculate_averaged_hamiltonian(x_z(i,:).',lambda_z(i,:).',params_zero);
end

fprintf('\n========================================\n');
fprintf('TEST 1: ZERO-THRUST TWO-BODY SANITY\n');
fprintf('========================================\n');
fprintf('max |p-p0| = %.15e\n', max(abs(x_z(:,1)-x0(1))));
fprintf('max |f-f0| = %.15e\n', max(abs(x_z(:,2)-x0(2))));
fprintf('max |g-g0| = %.15e\n', max(abs(x_z(:,3)-x0(3))));
fprintf('max |h-h0| = %.15e\n', max(abs(x_z(:,4)-x0(4))));
fprintf('max |k-k0| = %.15e\n', max(abs(x_z(:,5)-x0(5))));
fprintf('max |m-m0| = %.15e\n', max(abs(x_z(:,9)-x0(9))));
fprintf('max |H-H0| = %.15e\n', max(abs(H_z-H_z(1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2: Hamiltonian dynamics vs direct averaged state dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ydot_H = averaged_dynamics([x0; lambda0_paper_avg],params);
xdot_H = ydot_H(1:9);

xdot_direct = direct_average_xdot(x0,lambda0_paper_avg,params);

fprintf('\n========================================\n');
fprintf('TEST 2: HAMILTONIAN XDOT VS DIRECT AVG XDOT\n');
fprintf('========================================\n');
fprintf('norm(xdot_H - xdot_direct) = %.15e\n', norm(xdot_H - xdot_direct));
fprintf('relative error = %.15e\n', norm(xdot_H - xdot_direct)/max(1,norm(xdot_direct)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3: fixed paper-costate propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau,Y] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                [0 1],[x0; lambda0_paper_avg],odeopts);

x = Y(:,1:9);
lambda = Y(:,10:18);

H = zeros(length(tau),1);
sigma = zeros(length(tau),1);
T = zeros(length(tau),1);

for i = 1:length(tau)
    H(i) = calculate_averaged_hamiltonian(x(i,:).',lambda(i,:).',params);
    [sigma(i),T(i)] = thrust_history_simple(x(i,:).',lambda(i,:).',params);
end

res_fixed = terminal_residual(x(end,:).',lambda(end,:).',target);

fprintf('\n========================================\n');
fprintf('TEST 3: FIXED PAPER-COSTATE PROPAGATION\n');
fprintf('========================================\n');
fprintf('terminal residual norm = %.15e\n', norm(res_fixed));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 4: lsqnonlin correction from paper costates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lsqopts = optimoptions(@lsqnonlin, ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxIterations',100, ...
    'MaxFunctionEvaluations',2500);

lambda0_opt = lsqnonlin(@(lambda0) shooting_residual(lambda0,x0,target,params,odeopts), ...
                        lambda0_paper_avg, [], [], lsqopts);

[tau_opt,Y_opt] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                        [0 1],[x0; lambda0_opt],odeopts);

x_opt = Y_opt(:,1:9);
lambda_opt = Y_opt(:,10:18);

H_opt = zeros(length(tau_opt),1);
sigma_opt = zeros(length(tau_opt),1);
T_opt = zeros(length(tau_opt),1);

for i = 1:length(tau_opt)
    H_opt(i) = calculate_averaged_hamiltonian(x_opt(i,:).',lambda_opt(i,:).',params);
    [sigma_opt(i),T_opt(i)] = thrust_history_simple(x_opt(i,:).',lambda_opt(i,:).',params);
end

res_opt = terminal_residual(x_opt(end,:).',lambda_opt(end,:).',target);

fprintf('\n========================================\n');
fprintf('TEST 4: LSQNONLIN CORRECTED COSTATES\n');
fprintf('========================================\n');
fprintf('Optimized costates:\n');
disp(lambda0_opt);
fprintf('terminal residual norm = %.15e\n', norm(res_opt));
fprintf('Final p error = %.15e\n', x_opt(end,1)-target.p);
fprintf('Final f error = %.15e\n', x_opt(end,2)-target.f);
fprintf('Final g error = %.15e\n', x_opt(end,3)-target.g);
fprintf('Final h error = %.15e\n', x_opt(end,4)-target.h);
fprintf('Final k error = %.15e\n', x_opt(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda_opt(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda_opt(end,9));
fprintf('Final mass = %.15f\n', x_opt(end,9));
fprintf('max |H-H0| = %.15e\n', max(abs(H_opt-H_opt(1))));
fprintf('min sigma = %.15f\n', min(sigma_opt));
fprintf('max sigma = %.15f\n', max(sigma_opt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartesian plot of mean trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_cart = zeros(length(tau_opt),3);

for i = 1:length(tau_opt)
    s_mee = x_opt(i,1:6).';
    s_kep = MEE2KEP(s_mee);
    s_cart = KEP2CART(s_kep,params.mu);
    r_cart(i,:) = s_cart(1:3).';
end

lw = 1.5;

figure;
plot(tau_opt,x_opt(:,1),'LineWidth',lw); hold on;
plot(tau_opt,x_opt(:,2),'LineWidth',lw);
plot(tau_opt,x_opt(:,3),'LineWidth',lw);
plot(tau_opt,x_opt(:,4),'LineWidth',lw);
plot(tau_opt,x_opt(:,5),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('MEE states');
title('Corrected GTO-to-GEO Flat MEE State History');
legend('p','f','g','h','k','Location','best');

figure;
plot(tau_opt,lambda_opt(:,1),'LineWidth',lw); hold on;
plot(tau_opt,lambda_opt(:,2),'LineWidth',lw);
plot(tau_opt,lambda_opt(:,3),'LineWidth',lw);
plot(tau_opt,lambda_opt(:,4),'LineWidth',lw);
plot(tau_opt,lambda_opt(:,5),'LineWidth',lw);
plot(tau_opt,lambda_opt(:,6),'LineWidth',lw);
plot(tau_opt,lambda_opt(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('Costates');
title('Corrected GTO-to-GEO Flat Costate History');
legend('\lambda_p','\lambda_f','\lambda_g','\lambda_h','\lambda_k','\lambda_L','\lambda_m','Location','best');

figure;
plot(tau_opt,H_opt-H_opt(1),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('\tilde H-\tilde H(0)');
title('Averaged Hamiltonian Error');

figure;
plot(tau_opt,sigma_opt,'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('\sigma');
title('Thrust Modulation');

figure;
plot(tau_opt,T_opt,'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('T');
title('Thrust Magnitude');

figure;
plot(tau_opt,x_opt(:,9),'LineWidth',lw);
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
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = shooting_residual(lambda0,x0,target,params,odeopts)

    [~,Y] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                  [0 1],[x0; lambda0],odeopts);

    yf = Y(end,:).';
    xf = yf(1:9);
    lambdaf = yf(10:18);

    res = terminal_residual(xf,lambdaf,target);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test3.m
% Brute-force averaged shooting optimization with alpha optimization
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

params.q = 40;
params.epsilon_S = 1e-2;
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

target.p = p0;
target.f = f0;
target.g = g0;
target.h = 0;
target.k = 0;

lambda0_guess = [
  -0.464766811473970
  -9.147446698788391
   0.946304117863429
   7.231076914575487
  -1.437669479858417
  -0.013600630673835
   0.000000000006312
   0
   5.348790384121387
];

z0_guess = [lambda0_guess; alpha0];

odeopts = odeset('RelTol',1e-8,'AbsTol',1e-10);

optopts = optimoptions(@fminunc, ...
    'Display','iter', ...
    'Algorithm','quasi-newton', ...
    'OptimalityTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',10000);

z_opt = fminunc(@(z) shooting_cost_bruteforce_alpha(z,x0,target,params,odeopts), ...
                z0_guess,optopts);

lambda0_opt = z_opt(1:9);
alpha0_opt = z_opt(10);

x0(8) = alpha0_opt;

[tau,Y] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                [0 1],[x0; lambda0_opt],odeopts);

x = Y(:,1:9);
lambda = Y(:,10:18);

H = zeros(length(tau),1);
sigma = zeros(length(tau),1);
T = zeros(length(tau),1);

for i = 1:length(tau)
    H(i) = calculate_averaged_hamiltonian(x(i,:).',lambda(i,:).',params);
    [sigma(i),T(i)] = thrust_history_simple(x(i,:).',lambda(i,:).',params);
end

res = terminal_residual_alpha(x(end,:).',lambda(end,:).',target);

fprintf('\nOptimized initial costates:\n');
disp(lambda0_opt);

fprintf('alpha0_opt = %.15f\n', alpha0_opt);
fprintf('time of flight days = %.15f\n', alpha0_opt*TU/86400);
fprintf('terminal residual norm = %.15e\n', norm(res));
fprintf('Final p error = %.15e\n', x(end,1)-target.p);
fprintf('Final f error = %.15e\n', x(end,2)-target.f);
fprintf('Final g error = %.15e\n', x(end,3)-target.g);
fprintf('Final h error = %.15e\n', x(end,4)-target.h);
fprintf('Final k error = %.15e\n', x(end,5)-target.k);
fprintf('lambda_L(tf) = %.15e\n', lambda(end,6));
fprintf('lambda_m(tf) = %.15e\n', lambda(end,9));
fprintf('lambda_alpha(tf) = %.15e\n', lambda(end,8));
fprintf('Final mass = %.15f\n', x(end,9));
fprintf('max abs(H-H0) = %.15e\n', max(abs(H-H(1))));

r_cart = zeros(length(tau),3);

for i = 1:length(tau)
    s_mee = x(i,1:6).';
    s_kep = MEE2KEP(s_mee);
    s_cart = KEP2CART(s_kep,params.mu);
    r_cart(i,:) = s_cart(1:3).';
end

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
title('Optimized Brute-Force Averaged MEE States');
legend('p','f','g','h','k','Location','best');

figure;
plot(tau,lambda(:,1),'LineWidth',lw); hold on;
plot(tau,lambda(:,2),'LineWidth',lw);
plot(tau,lambda(:,3),'LineWidth',lw);
plot(tau,lambda(:,4),'LineWidth',lw);
plot(tau,lambda(:,5),'LineWidth',lw);
plot(tau,lambda(:,6),'LineWidth',lw);
plot(tau,lambda(:,8),'LineWidth',lw);
plot(tau,lambda(:,9),'LineWidth',lw);
grid on; box on;
xlabel('\tau');
ylabel('Costates');
title('Optimized Brute-Force Averaged Costates');
legend('\lambda_p','\lambda_f','\lambda_g','\lambda_h','\lambda_k','\lambda_L','\lambda_\alpha','\lambda_m','Location','best');

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

function J = shooting_cost_bruteforce_alpha(z,x0,target,params,odeopts)

    lambda0 = z(1:9);
    alpha0 = z(10);

    if alpha0 <= 0
        J = 1e30;
        return
    end

    x0(8) = alpha0;

    [~,Y] = ode45(@(tau,y) averaged_dynamics(y,params), ...
                  [0 1],[x0; lambda0],odeopts);

    yf = Y(end,:).';
    xf = yf(1:9);
    lambdaf = yf(10:18);

    res = terminal_residual_alpha(xf,lambdaf,target);

    W = diag([1 1 1 10 10 1 1 1]);

    J = res.'*W*res;

end

function res = terminal_residual_alpha(xf,lambdaf,target)

    res = [
        xf(1) - target.p
        xf(2) - target.f
        xf(3) - target.g
        xf(4) - target.h
        xf(5) - target.k
        lambdaf(6)
        lambdaf(9)
        lambdaf(8)
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
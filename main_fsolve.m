%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 568 Final Project
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc; clear; close all;
% 
% u_max = 0.1;
%%
%% BVP
% Constants

mu = 1; % Sun's gravitational parameter (km^3/s^2)
l = 1; % [1 a.u]
t = 1;
v = 1;
a = 1;
tf = 8;
u_max = 0.1;
m_nd = 1; % 1500kg
Isp = 7.9639e-4; %s
g0 =  1.6543e3; 

% Nondimensionalizing constants

% Initial state [r0, v0] (Example: circular orbit at 7000 km altitude)
a_earth_0 = 1; %[A.U.]
nu_earth_0 = 0; %[rad]

a_mars_0 = 1.524; %[A.U.]
nu_mars_0 = pi; %[rad]

r_earth_0 = [a_earth_0; 0];
r0 = r_earth_0;
v0 = [0; sqrt(mu / norm(r_earth_0))]; 

% Initial guess for costates
lambda0_guess = [0.849671262182387; 0.101610652998729; 0.201818072190866; 1.084649182385908; 0.5]; % This is for rho = 1, initial
rho = 1;

options =  optimoptions(@fsolve,'OptimalityTolerance', 1e-12, 'FunctionTolerance', 1e-12); % CHANGE AGAIN LATER

% Solve for optimal costates using fsolve

lambda0_opt = fsolve(@(lambda0) shooting_function_averaged(lambda0, r0, v0, m_nd, mu, tf, r_mars_0, v0_mars, u_max, rho, Isp, g0), lambda0_guess, options);

%Also need to solve for tf with this
tspan = [0 tf];


% Integrate with optimal costates
options2 = odeset("RelTol",1e-12, "AbsTol",1e-12);
[t, Y] = ode45(@(t, y) indirect_orbit_dynamics_mid(y, mu, u_max, rho, Isp, g0), tspan, [r0; v0; m_nd; lambda0_opt], options2);

r = Y(:,1:2);
v = Y(:,3:4);
m = Y(:, 5);
x = [r, v, m];

lambda = Y(:,6:10);
% Calculating control history
u = optimal_controls(m, lambda, u_max, rho, Isp, g0)';

% % Calculating Hamiltonian
% r_norm = sqrt(r(:,1).^2 + r(:,2).^2);
H = zeros(numel(t));
for i = 1:numel(t)
    H(i) = calculate_hamiltonian(x(i, :)', lambda(i, :)', mu, u_max, rho, Isp, g0);
end

H0 = H(1);

u_norms = vecnorm(u);
fuel = sum(u_norms);
%% Calculating mars orbit
[t_mars, X_mars] = ode45(@(t, X_mars) orbit_dynamics(t, X_mars, mu), tspan, [r_mars_0; v0_mars]);
[t_earth, X_earth] = ode45(@(t, X_earth) orbit_dynamics(t, X_earth, mu), tspan, [r_earth_0; v0]);
% Maybe it's better to just do a circular orbit here


%% Plotting
% Trajectories
figure
plot(Y(:,1), Y(:,2))
xlabel("X (A.U)")
ylabel("Y (A.U)")
hold on
plot(X_mars(:,1), X_mars(:,2))
plot(X_earth(:, 1), X_earth(:, 2))
legend({"Rocket Trajectory", "Mars Trajectory", "Earth Trajectory"}, 'Location', 'northwest')
grid on
axis equal
title("Trajectory")

% Costate History
figure
plot(t, lambda(:, 1))
hold on
plot(t, lambda(:, 2))
plot(t, lambda(:, 3))
plot(t, lambda(:, 4))
plot(t, lambda(:, 5))
grid on
xlabel("time (n.d)")
ylabel("Costates")
legend("\lambda_rx", "\lambda_ry", "\lambda_vx", "\lambda_vy", "\lambda_m")
title("Costate History")

% Control History
figure
plot(t, u(1,:))
hold on
plot(t, u(2,:))
ylabel("Control u(t)")
xlabel("Time (n.d.)")
grid on
title("Control History")

figure
plot(t, sqrt(u(1,:).^2 + u(2,:).^2))
ylabel("Norm of u(t)")
xlabel("Time (n.d.)")
grid on
title("Control Norm")

disp(lambda0_opt)
% Hamiltonian plot
figure
plot(t, H(:,1)-H0, "r");
xlabel("Time (n.d)")
ylabel("H - H0")
grid on
title("Hamiltonian History")

figure
plot(t, m)
title("Mass History")
xlabel("Time (n.d.)")
ylabel("Mass (n.d.)")
grid on

%% User Defined Function
function F = shooting_function_averaged(lambda0_guess, x0_mean, params)
    tauspan = [0 1];
    [~, Y] = ode45(@(tau, y) averaged_dynamics(y, params), tauspan, [x0_mean; lambda0_guess], options);
    
    yf = Y(end, :)';
    xf = yf(1:9);
    lambdaf = yf(10:18);
    
    p_f = xf(1); f_f = xf(2); g_f = xf(3); h_f = xf(4); k_f = xf(5);
    
    p_target = 42165; 
    ecc_target = 0;
    inc_target = 0; % Radians
    
    % 1-5: The Orbital Elements (Fully Constrained)
    F(1) = p_f - p_target; % Semilatus
    F(2) = f_f - f_target; 
    F(3) = g_f - g_target;
    F(4) = h_f - h_target;
    F(5) = k_f - k_target;
    
    % We don't care where the spacecraft is in the orbit when the transfer ends.
    F(6) = lambdaf(6); 
    F(7) = lambdaf(9); 
    F(8) = lambdaf(7); 
    F(9) = lambdaf(8);
end


function [u] = optimal_controls(m, lambda, u_max, rho, Isp, g0)
    lambda_v = lambda(:, 3:4);
    lambda_m = lambda(:, 5);

    p = - lambda_v ./ m;
    p_norm = vecnorm(p');
    S = p_norm' - (1 - lambda_m / (Isp * g0));
    Gamma = (u_max / 2) .* (1 + tanh(S / rho));
    u = Gamma .* p ./ p_norm'; % p / p_norm is u* hat, with Gamma being magnitude
end
function [u] = optimal_control(m, lambda, u_max, rho, Isp, g0)
    lambda_v = lambda(3:4);
    lambda_m = lambda(5);

    p = - lambda_v ./ m;
    p_norm = vecnorm(p');
    S = p_norm' - (1 - lambda_m / (Isp * g0));
    Gamma = (u_max / 2) .* (1 + tanh(S / rho));
    u = Gamma .* p ./ p_norm'; % p / p_norm is u* hat, with Gamma being magnitude
end
function [H] = calculate_hamiltonian(x, lambda, mu, u_max, rho, Isp, g0)
    m = x(5);
    u = optimal_control(m, lambda, u_max, rho, Isp, g0)';
    dxdt = indirect_orbit_dynamics_mid([x; lambda], mu, u_max, rho, Isp, g0);
    H = vecnorm(u) + lambda.' * dxdt(1:5);
end
function [xdot] = state_equation(x, lambda, mu, B, u_max, rho)
    rvec = x(1:2);
    vvec = x(3:4);
    r = norm(rvec);
    
    f0 = [vvec; -mu / r ^ 3 * rvec];

    u = optimal_control(lambda, B, u_max, rho);
    xdot = f0 + B * u;
end
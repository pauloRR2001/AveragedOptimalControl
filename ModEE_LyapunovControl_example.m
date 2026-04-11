clc; clear; close all;

%% Constants and scaling
mu_km = 398600.4418;             % [km^3/s^2]

%% Initial Keplerian elements of deputy (spacecraft)
a_d0     = 10000;   % km
ecc_d0     = 0.3;
inc_d0     = deg2rad(10);
raan_d0 = deg2rad(0);
argp_d0 = deg2rad(0);
M_d0     = 0;       % true anomaly
nu_d0 = M_d0;

%% Keplerian elements of target
a_t     = 60000;  % km
ecc_t     = 0.7;
inc_t     = deg2rad(130);
raan_t = deg2rad(180);
argp_t = deg2rad(270);

% Target orbital elements (fixed)
[p_t,f_t,g_t,hf_t,k_t,L_t] = kep2modEqui(a_t,ecc_t,inc_t,raan_t,argp_t,0);

L_star = p_t;                  % [km]
mu = 1;                          % nondimensional
t_star = sqrt(L_star^3 / mu_km); % [s]

[p0,f0,g0,hf0,k0,L0] = kep2modEqui(a_d0,ecc_d0,inc_d0,raan_d0,argp_d0,nu_d0);

% Nondimensionalize initial state
x0_d = [p0; f0; g0; hf0; k0; L0];
x0_d(1) = x0_d(1) / L_star;  % nondimensionalize semi-major axis

x_t_slow = [p_t; f_t; g_t; hf_t; k_t];
x_t_slow(1) = x_t_slow(1) / L_star;

%% Time span
%T_c_dim = 2*pi*sqrt(a_d0^3 / mu_km); % [s] chief period (based on initial orbit)
%T_c = T_c_dim / t_star;           % nondimensional
day = 1*3600*24 / t_star;
tspan = [0, 10*day];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Control design
K = eye(5);              % Lyapunov gain
P = 1 * eye(5);          % PD decay rate
u_max = 0.1/1000 * t_star^2 / L_star*10^10; % nondim bound for .1 m/s^2

%% Integrate
[t, x_d] = ode45(@(t, x_d) ...
    controlled_orbit_dynamics_equinoctial(t, x_d, x_t_slow, K, P, u_max, mu), ...
    tspan, x0_d, opts);

% Optional: convert back to Cartesian for plotting
% x_d_kep = elementArrayConvert([x_d(:,1)*L_star, x_d(:,2:6)],'MOD','KEP',mu_km);
% x_d_cart = elementArrayConvert(x_d_kep,'KEP','CAR',mu_km);
for j=1:length(t)
    x_d_cart(j,1:6) = MEE2cart([x_d(j,1)*L_star, x_d(j,2:6)],mu_km);

end
x_d_kep = elementArrayConvert(x_d_cart,'CAR','KEP',mu_km);

% Plot results
figure;
plot3(x_d_cart(:,1), x_d_cart(:,2), x_d_cart(:,3), 'k','LineWidth', 1.0);
xlabel('x [km]','interpreter','latex'); ylabel('y [km]','interpreter','latex');
zlabel('z [km]','interpreter','latex');
title('Deputy Trajectory (MEE Control)','interpreter','latex');
grid on; axis equal;
hold on
earthy(6378.137,'Earth',1,[0 0 0])
plotOrbit3(raan_t, inc_t, argp_t, a_t*(1-ecc_t^2), ecc_t, linspace(0,2*pi,1000), 'r', 1, 1, [0,0,0],0,1.2)
plot3(x_d_cart(1,1), x_d_cart(1,2), x_d_cart(1,3), 'bo','LineWidth', 1.2);
plot3(x_d_cart(end,1), x_d_cart(end,2), x_d_cart(end,3), 'rx','LineWidth', 1.2);
legend('Deputy Trajectory','','Target Orbit','IC','FC','interpreter','latex')

%% plots

for j=1:length(t)
    del_x(j,1:6) = x_d(j,:)-[p_t/L_star,f_t,g_t,hf_t,k_t,L_t];
end
plot_relative_keplerian_elements(t*t_star, [del_x(:,1)*L_star, del_x(:,2:end)])

% controller recreation
% x0_d = [a_d0; ecc_d0; inc_d0; raan_d0; argp_d0; M_d0];
for j=1:length(t)
    B = gauss_variational_matrix(x_d(j,1), x_d(j,2), x_d(j,3), x_d(j,5), x_d(j,4)...
        , x_d(j,6), mu);
    
    % Take only first 5 rows (slow variables)
    B_slow = B(1:5, :);  % [5 x 3]
    
    % Compute control law
    y = B_slow' * K' * P * del_x(j,1:5)';
    u(j,1:3) = - (B_slow' * K' * K * B_slow) \ y;

    if norm(u(j,1:3)) > u_max
        u(j,1:3) = u_max * u(j,1:3) / norm(u(j,1:3));
    end
    u_norm(j) = norm(u(j,1:3));
end
DV = trapz(t,u_norm);
disp(DV*L_star/t_star)

figure
plot(t*t_star/(3600*24),u*L_star/t_star^2,'LineWidth',1.2)
hold on
plot(t*t_star/(3600*24),u_norm*L_star/t_star^2,'k','LineWidth',1.2)
%yline(u_max*L_star/t_star^2,'r')
%yline(-u_max*L_star/t_star^2,'r')
grid on
title('Evolution of Control Acceleration Input U','interpreter','latex')
xlabel('Time [days]','interpreter','latex')
ylabel('u $[km/s^2]$','interpreter','latex')
legend('$u_r$','$u_{\theta}$','$u_h$','$||u||_2$','$u_{max}$','','interpreter','latex')

%% func

% Dynamics function
function xdot = controlled_orbit_dynamics_equinoctial(~, x_d, x_t_slow, K, P, u_max, mu)
    % Unpack MEE state
    p_d  = x_d(1);
    f_d  = x_d(2);
    g_d  = x_d(3);
    h_d  = x_d(4);
    k_d  = x_d(5);
    L_d  = x_d(6);

    % Compute error in slow states
    x_slow_d = x_d(1:5);
    del_x_slow = x_slow_d - x_t_slow;

    % Gauss variational matrix in MEE
    B = gauss_variational_matrix(p_d, f_d, g_d, h_d, k_d, L_d, mu);
    B_slow = B(1:5, :);  % only use rows for p,f,g,h,k

    % Lyapunov PD-like control law
    y = B_slow' * K' * P * del_x_slow;
    u = - (B_slow' * K' * K * B_slow) \ y;

    % Saturate control input
    if norm(u) > u_max
        u = u_max * u / norm(u);
    end

    % Dynamics update
    xdot = B * u;

    % Add mean longitude rate (natural drift)
    q_d = 1 + f_d * cos(L_d) + g_d * sin(L_d);
    n_d = sqrt(mu / p_d^3) * q_d^2;
    xdot(6) = xdot(6) + n_d;
end


function B = gauss_variational_matrix(p, f, g, h, k, L, mu)
% GAUSS_VARIATIONAL_MATRIX  Computes the Gauss matrix B(x) for mapping
% accelerations in the RTN frame to rates of change of Keplerian elements.
%
% Inputs:
%   a      - semi-major axis (nondimensional)
%   e      - eccentricity
%   i      - inclination [rad]
%   omega  - argument of periapsis [rad]
%   Omega  - RAAN [rad]
%   M      - mean anomaly [rad]
%   mu     - gravitational parameter (nondimensional)
%
% Output:
%   B      - 6x3 Gauss variational matrix

    q = 1+f*cos(L)+g*sin(L);

    s = sqrt(1+h^2+k^2);

    B = sqrt(p/mu)*[0, 2*p/q, 0;...
        sin(L), ((q+1)*cos(L)+f)/q, -g*(h*sin(L)-k*cos(L))/q;...
        -cos(L), ((q+1)*sin(L)+g)/q, f*(h*sin(L)-k*cos(L))/q;...
        0, 0, (s^2/(2*q))*cos(L);...
        0, 0, (s^2/(2*q))*sin(L);...
        0, 0, (h*sin(L)-k*cos(L))/q];
end

function plot_relative_keplerian_elements(t, delta_x_kep)
% Plots relative Keplerian elements over time using LaTeX formatting.
%
% Inputs:
%   t              : time vector [s]
%   delta_x_kep    : matrix of size (N x 6) containing:
%                    [δa (km), δe, δi (rad), δΩ (rad), δω (rad), δM (rad)]

    % Extract relative elements
    %[p_t,f_t,g_t,hf_t,k_t,L_t]
    dp    = delta_x_kep(:,1);   % [km]
    df    = delta_x_kep(:,2);   % [-]
    dg    = (delta_x_kep(:,3));   % [deg]
    dh = (delta_x_kep(:,4));   % [deg]
    dk = (delta_x_kep(:,5));   % [deg]
    dL    = rad2deg(delta_x_kep(:,6));   % [deg]

    % Plotting
    figure;
    subplot(3,2,1);
    plot(t/(3600*24), dp, 'k', 'LineWidth', 1.3);
    xlabel('Time [days]', 'Interpreter', 'latex');
    ylabel('$\delta p$ [km]', 'Interpreter', 'latex');
    title('Relative Semi-Latus Rectum', 'Interpreter', 'latex');
    hold on
    yline(0,'r')
    grid on;

    subplot(3,2,2);
    plot(t/(3600*24), df, 'k', 'LineWidth', 1.3);
    xlabel('Time [days]', 'Interpreter', 'latex');
    ylabel('$\delta f$', 'Interpreter', 'latex');
    title('Relative f Element', 'Interpreter', 'latex');
    grid on;
    hold on
    yline(0,'r')

    subplot(3,2,3);
    plot(t/(3600*24), dg, 'k', 'LineWidth', 1.3);
    xlabel('Time [days]', 'Interpreter', 'latex');
    ylabel('$\delta g$', 'Interpreter', 'latex');
    title('Relative g Element', 'Interpreter', 'latex');
    grid on;
    hold on
    yline(0,'r')

    subplot(3,2,4);
    plot(t/(3600*24), dh, 'k', 'LineWidth', 1.3);
    xlabel('Time [days]', 'Interpreter', 'latex');
    ylabel('$\delta h$', 'Interpreter', 'latex');
    title('Relative h Element', 'Interpreter', 'latex');
    grid on;
    hold on
    yline(0,'r')

    subplot(3,2,5);
    plot(t/(3600*24), dk, 'k', 'LineWidth', 1.3);
    xlabel('Time [days]', 'Interpreter', 'latex');
    ylabel('$\delta k$', 'Interpreter', 'latex');
    title('Relative k Element', 'Interpreter', 'latex');
    grid on;
    hold on
    yline(0,'r')

    subplot(3,2,6);
    plot(t/(3600*24), dL, 'k', 'LineWidth', 1.3);
    xlabel('Time [days]', 'Interpreter', 'latex');
    ylabel('$\delta L$ [deg]', 'Interpreter', 'latex');
    title('Relative True Longitude', 'Interpreter', 'latex');
    grid on;

    sgtitle('Time Evolution of Relative Equinoctial Elements', 'Interpreter', 'latex');
end

function x_cart = MEE2cart(x_MEE, mu)
    % Extract elements
    p = x_MEE(1);
    f = x_MEE(2);
    g = x_MEE(3);
    h = x_MEE(4);
    k = x_MEE(5);
    L = x_MEE(6);

    % Compute auxiliary quantities
    s2 = 1 + h^2 + k^2;
    r = p / (1 + f*cos(L) + g*sin(L));
    root_mu_p = sqrt(mu / p);

    % Position vector
    x = r/s2 * (cos(L) + (s2 - 1)*cos(L) + 2*h*k*sin(L));
    y = r/s2 * (sin(L) - (s2 - 1)*sin(L) + 2*h*k*cos(L));
    z = r/s2 * (2*(h*sin(L) - k*cos(L)));

    % Velocity vector
    vx = -1/s2 * root_mu_p * ( ...
        sin(L) + (s2 - 1)*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + (s2 - 1)*g );

    vy = -1/s2 * root_mu_p * ( ...
       -cos(L) + (s2 - 1)*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + (s2 - 1)*f );

    vz = -1/s2 * root_mu_p * ( ...
        -2*(h*cos(L) - k*sin(L)) + f*h + g*k );

    x_cart = [x, y, z, vx, vy, vz];
end

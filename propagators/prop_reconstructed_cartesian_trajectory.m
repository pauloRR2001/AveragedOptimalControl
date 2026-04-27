function out = prop_reconstructed_cartesian_trajectory(tau_avg, x_avg, lambda_avg, params, npts_per_orbit)
% Reconstructs a realistic-looking controlled trajectory by propagating the
% UNAVERAGED MEE dynamics densely, while using the converged AVERAGED
% costates as a feedback law for thrust direction and throttle.
%
% INPUTS
%   tau_avg        [N x 1] averaged independent-variable grid
%   x_avg          [N x 9] averaged state history
%   lambda_avg     [N x 9] averaged costate history
%   params         struct with fields used by xdot_unavg_sym_fun
%   npts_per_orbit optional density target for reconstruction (default 200)
%
% REQUIRED FUNCTIONS ALREADY IN YOUR PATH
%   xdot_unavg_sym_fun
%   get_MEE_B_matrix
%   MEE2KEP
%   KEP2CART
%
% OUTPUT
%   out struct with dense reconstructed histories:
%       .tau
%       .x
%       .lambda
%       .sigma
%       .T
%       .u_rth
%       .r_cart
%       .v_cart
%
% USAGE
%   out = plot_reconstructed_cartesian_trajectory(tau, x, lambda, params);
%
% IDEA
%   The averaged solution only captures slow secular evolution, so plotting
%   its sparse points in Cartesian space gives straight chords. This
%   function restores the fast orbital motion by propagating the unaveraged
%   state dynamics densely, while using interpolated averaged costates to
%   generate the PMP thrust law at each dense step.

    if nargin < 5 || isempty(npts_per_orbit)
        npts_per_orbit = 200;
    end

    tau_avg = tau_avg(:);
    if size(x_avg,1) ~= numel(tau_avg) || size(lambda_avg,1) ~= numel(tau_avg)
        error('tau_avg, x_avg, and lambda_avg must have consistent lengths.');
    end
    if size(x_avg,2) ~= 9 || size(lambda_avg,2) ~= 9
        error('x_avg and lambda_avg must be N-by-9.');
    end

    tau0 = tau_avg(1);
    tauf = tau_avg(end);
    x0 = x_avg(1,:).';

    % Estimate number of reconstructed points from longitude accumulation.
    dL = diff(unwrap(x_avg(:,6)));
    nrev_est = max(1, ceil(abs(sum(dL))/(2*pi)));
    ndense = max(numel(tau_avg), nrev_est*npts_per_orbit);

    tau_dense_eval = linspace(tau0, tauf, ndense).';

    odeopts = odeset('RelTol',1e-10,'AbsTol',1e-11);

    % Dense unaveraged propagation using interpolated averaged costates.
    [tau_dense, x_dense] = ode113(@(tau,x) unavg_rhs_with_interp_costate(tau, x, tau_avg, lambda_avg, params), ...
                                  tau_dense_eval, x0, odeopts);

    % Interpolate costates on the dense grid for diagnostics / output.
    lambda_dense = interp_costates(tau_dense, tau_avg, lambda_avg);

    % Recover thrust quantities.
    sigma = zeros(numel(tau_dense),1);
    T = zeros(numel(tau_dense),1);
    u_rth = zeros(numel(tau_dense),3);

    for i = 1:numel(tau_dense)
        [sigma(i), T(i), u_rth(i,:)] = control_from_costates(x_dense(i,:).', lambda_dense(i,:).', params);
    end

    % Convert dense MEE history to Cartesian.
    r_cart = zeros(numel(tau_dense),3);
    v_cart = zeros(numel(tau_dense),3);

    for i = 1:numel(tau_dense)
        s_mee = x_dense(i,1:6).';
        s_kep = MEE2KEP(s_mee);
        s_cart = KEP2CART(s_kep, params.mu);
        r_cart(i,:) = s_cart(1:3).';
        v_cart(i,:) = s_cart(4:6).';
    end

    out.tau = tau_dense;
    out.x = x_dense;
    out.lambda = lambda_dense;
    out.sigma = sigma;
    out.T = T;
    out.u_rth = u_rth;
    out.r_cart = r_cart;
    out.v_cart = v_cart;

end

function xdot = unavg_rhs_with_interp_costate(tau, x, tau_avg, lambda_avg, params)

    lambda = interp_costates(tau, tau_avg, lambda_avg).';
    Lq = x(6);

    xdot = xdot_unavg_sym_fun( ...
        x, lambda, Lq, ...
        params.mu, params.c, params.T_min, params.T_max, ...
        params.eps_S, params.RE, params.J2);

end

function lambda_q = interp_costates(tau_q, tau_avg, lambda_avg)

    tau_q = tau_q(:);
    lambda_q = zeros(numel(tau_q), 9);

    for j = 1:9
        lambda_q(:,j) = interp1(tau_avg, lambda_avg(:,j), tau_q, 'pchip', 'extrap');
    end

end

function [sigma, T, u_rth] = control_from_costates(x, lambda, params)

    B = get_MEE_B_matrix(x, params.mu);

    lam_MEE = lambda(1:6);
    lam_m = lambda(9);

    pvec = B.' * lam_MEE;
    pnorm = norm(pvec);

    if pnorm > 1e-14
        u_rth = (-pvec / pnorm).';
    else
        u_rth = [0 0 0];
    end

    S = -(params.c / x(9)) * pnorm - lam_m + 1;
    sigma = 0.5 * (1 - S / sqrt(S^2 + params.eps_S^2));
    T = params.T_min + (params.T_max - params.T_min) * sigma;

end
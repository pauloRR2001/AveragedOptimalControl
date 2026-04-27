function out = prop_kepler_orbit_from_mee_for_time(s_mee0, t_f, mu, npts)
% Plots the two-body orbit corresponding to an initial MEE state by:
%   1) converting MEE -> Keplerian once at t = 0
%   2) advancing mean anomaly as M(t) = M0 + n t
%   3) converting each Keplerian sample to Cartesian
%
% INPUTS
%   s_mee0  [6x1] or [1x6] = [p f g h k L] at t = 0
%   t_final scalar duration in the same time units as mu
%   mu      gravitational parameter
%   npts    optional number of points, default 2000
%
% ASSUMPTION
%   Two-body Keplerian motion only. No thrust, no J2, no perturbations.
%
% OUTPUT
%   out struct with fields:
%       .t
%       .kep
%       .cart
%
% NOTE
%   This requires:
%       MEE2KEP.m
%       KEP2CART.m
%
% EXAMPLE
%   out = plot_kepler_orbit_from_mee_for_time([p f g h k L0], 5*86400/TU, 1, 3000);

    if nargin < 4 || isempty(npts)
        npts = 2000;
    end

    s_mee0 = s_mee0(:);
    if numel(s_mee0) ~= 6
        error('s_mee0 must be a 6-element MEE state [p f g h k L].');
    end

    s_kep0 = MEE2KEP(s_mee0);

    a    = s_kep0(1);
    ecc  = s_kep0(2);
    inc  = s_kep0(3);
    raan = s_kep0(4);
    argp = s_kep0(5);
    M0   = s_kep0(6);

    if a <= 0
        error('Semimajor axis must be positive.');
    end

    if ecc >= 1
        error('This function assumes an elliptic orbit: ecc must be < 1.');
    end

    n = sqrt(mu / a^3);
    tau = 2*pi/n;

    t = linspace(0, t_f, npts).';
    kep_hist = zeros(npts, 6);
    cart_hist = zeros(npts, 6);

    for i = 1:npts
        M = wrapTo2Pi_local(M0 + n*t(i));
        s_kep = [a; ecc; inc; raan; argp; M];
        s_cart = KEP2CART(s_kep, mu);

        kep_hist(i,:) = s_kep.';
        cart_hist(i,:) = s_cart.';
    end

    out.t = t;
    out.kep = kep_hist;
    out.cart = cart_hist;

end

function ang = wrapTo2Pi_local(ang)
    ang = mod(ang, 2*pi);
    if ang < 0
        ang = ang + 2*pi;
    end
end
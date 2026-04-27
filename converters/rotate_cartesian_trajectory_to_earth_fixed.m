function out = rotate_cartesian_trajectory_to_earth_fixed(t, r_eci, omega_E, theta0)
%ROTATE_CARTESIAN_TRAJECTORY_TO_EARTH_FIXED Rotate inertial positions into an Earth-rotating frame.
%
%   out = rotate_cartesian_trajectory_to_earth_fixed(t, r_eci, omega_E)
%   out = rotate_cartesian_trajectory_to_earth_fixed(t, r_eci, omega_E, theta0)
%
% INPUTS
%   t        [N x 1] or [1 x N]
%            Time history in nondimensional time units.
%
%   r_eci    [N x 3]
%            Cartesian position history in the inertial frame.
%            Only position is rotated. No velocity transport terms are used.
%
%   omega_E scalar
%            Earth rotation rate in nondimensional units, rad / TU.
%
%   theta0  scalar, optional
%            Initial Earth rotation angle at t(1), in radians.
%            Default is 0.
%
% OUTPUT
%   out struct with fields:
%       .t          [N x 1] time vector
%       .theta      [N x 1] Earth rotation angle history
%       .r_rot      [N x 3] position history in rotating frame
%       .r_eci      [N x 3] original inertial position history
%       .omega_E    scalar Earth rotation rate used
%       .theta0     scalar initial angle used
%
% CONVENTION
%   The rotating-frame position is computed as
%
%       r_rot = R3(theta)^T * r_eci
%
%   where R3(theta) rotates vectors from the rotating frame to the inertial
%   frame by angle theta about +z. Therefore R3(theta)^T maps inertial
%   vectors into the rotating frame.
%
%   Explicitly:
%
%       r_rot = [ cos(theta)   sin(theta)  0
%                -sin(theta)   cos(theta)  0
%                  0            0          1 ] * r_eci
%
%   No transport theorem terms are applied because this function only
%   rotates positions, not velocities.

    if nargin < 4 || isempty(theta0)
        theta0 = 0;
    end

    t = t(:);

    if size(r_eci,2) ~= 3
        error('r_eci must be an [N x 3] position array.');
    end

    if size(r_eci,1) ~= numel(t)
        error('Length of t must match the number of rows in r_eci.');
    end

    if ~isscalar(omega_E)
        error('omega_E must be a scalar in rad / nondimensional time unit.');
    end

    theta = theta0 + omega_E*(t - t(1));

    r_rot = zeros(size(r_eci));

    for i = 1:numel(t)

        c = cos(theta(i));
        s = sin(theta(i));

        C_I_to_R = [
             c,  s,  0
            -s,  c,  0
             0,  0,  1
        ];

        r_rot(i,:) = (C_I_to_R * r_eci(i,:).').';

    end

    out.t = t;
    out.theta = theta;
    out.r_rot = r_rot;
    out.r_eci = r_eci;
    out.omega_E = omega_E;
    out.theta0 = theta0;

end
% -------------------------------------------------------------------------
% Function: cart2kepECC
% Author: [Paulo Ramirez]
% Date: [2/12/2025]
%
% Description:
% Converts Cartesian state vectors (position and velocity) into orbital 
% elements using eccentricity vector components. This is useful for 
% near-circular orbits where the argument of perigee is ill-defined.
%
% Inputs:
%   x, y, z     - Cartesian position components [km]
%   vx, vy, vz  - Cartesian velocity components [km/s]
%   mu          - Standard gravitational parameter of Earth [km^3/s^2]
%
% Outputs:
%   a      - Semi-major axis [km]
%   eccx   - X-component of eccentricity (dimensionless)
%   eccy   - Y-component of eccentricity (dimensionless)
%   inc    - Inclination [radians]
%   raan   - Right Ascension of the Ascending Node (RAAN) [radians]
%   theta  - Argument of latitude [radians]
% -------------------------------------------------------------------------

function [a, eccx, eccy, inc, raan, theta] = cart2kepECC(x, y, z, vx, vy, vz, mu)

    % Construct position and velocity vectors
    r_vec_xyz = [x, y, z]; % Position vector (km)
    v_vec_xyz = [vx, vy, vz]; % Velocity vector (km/s)

    % Compute radial distance (norm of position vector)
    r_mag = norm(r_vec_xyz);
    
    % Compute unit radial vector
    r_hat_xyz = r_vec_xyz / r_mag;

    % Compute velocity magnitude
    v_mag = norm(v_vec_xyz);

    % Compute specific angular momentum vector
    h_vec_xyz = cross(r_vec_xyz, v_vec_xyz);
    h = norm(h_vec_xyz); % Angular momentum magnitude
    
    % Compute unit angular momentum vector
    h_hat_xyz = h_vec_xyz / h;

    % Compute unit true anomaly direction (perpendicular to h and r)
    theta_hat_xyz = cross(h_hat_xyz, r_hat_xyz);

    % Transformation matrix from **ECI to RTH frame**
    ICR = [r_hat_xyz', theta_hat_xyz', h_hat_xyz']; 
    RCI = ICR'; % Transpose to get RTH to ECI transformation

    % Transform velocity vector to RTH frame
    v_vec_rth = (RCI * v_vec_xyz')'; % Transpose to maintain row vector format

    % Compute specific orbital energy
    eps = 0.5 * v_mag^2 - mu / r_mag;

    % Compute semi-major axis
    a = -mu / (2 * eps);

    % Compute semi-latus rectum
    p = (h^2) / mu;

    % Compute eccentricity magnitude
    ecc = sqrt(1 - p / a);

    % Compute true anomaly based on radial velocity sign
    if v_vec_rth(1) >= 0
        f = acos((p / r_mag - 1) / ecc);
    else
        f = -acos((p / r_mag - 1) / ecc);
    end

    %% Compute 3D Orbital Elements

    % Compute inclination
    inc = acos(h_hat_xyz(3));

    % Compute Right Ascension of the Ascending Node (RAAN)
    raan_candidates = [asin(h_hat_xyz(1) / sin(inc)), ...
                      pi - asin(h_hat_xyz(1) / sin(inc)), ...
                      acos(-h_hat_xyz(2) / sin(inc)), ...
                      2*pi - acos(-h_hat_xyz(2) / sin(inc))];

    raan = nonuniqueAngle(raan_candidates); % Resolve angle ambiguity

    % Compute argument of latitude (theta)
    theta_candidates = [asin(r_hat_xyz(3) / sin(inc)), ...
                        pi - asin(r_hat_xyz(3) / sin(inc)), ...
                        acos(theta_hat_xyz(3) / sin(inc)), ...
                        2*pi - acos(theta_hat_xyz(3) / sin(inc))];

    theta = nonuniqueAngle(theta_candidates); % Resolve angle ambiguity

    % Compute argument of perigee
    argp = theta - f;

    % Compute eccentricity components
    eccx = ecc * cos(argp);
    eccy = ecc * sin(argp);
    
end

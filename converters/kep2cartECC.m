% -------------------------------------------------------------------------
% Function: kep2cartECC
% Author: [Paulo Ramirez]
% Date: [2/12/2025]
%
% Description:
% Converts orbital elements (using eccentricity components) into 
% Cartesian state vectors (position and velocity). This version accounts
% for near-circular orbits where the argument of perigee is ambiguous.
%
% Inputs:
%   a      - Semi-major axis [km]
%   eccx   - X-component of eccentricity (dimensionless)
%   eccy   - Y-component of eccentricity (dimensionless)
%   inc    - Inclination [radians]
%   raan   - Right Ascension of the Ascending Node (RAAN) [radians]
%   theta  - True anomaly [radians]
%   mu     - Gravitational parameter [km^3/s^2]
%
% Outputs:
%   x, y, z    - Position in Cartesian coordinates [km]
%   vx, vy, vz - Velocity in Cartesian coordinates [km/s]
% -------------------------------------------------------------------------

function [x, y, z, vx, vy, vz] = kep2cartECC(a, eccx, eccy, inc, raan, theta, mu)
    
    % Compute total eccentricity magnitude
    ecc = sqrt(eccx^2 + eccy^2);

    % Compute argument of perigee (argp) considering quadrant checks
    % Possible values from both cosine and sine definitions
    argp = [acos(eccx / ecc), 2*pi - acos(eccx / ecc), asin(eccy / ecc), pi - asin(eccy / ecc)];
    
    % Resolve ambiguity in argument of perigee selection
    argp = nonuniqueAngle(argp);

    % Compute semi-latus rectum (p) and angular momentum (h)
    p = a * (1 - ecc^2);
    h = sqrt(p * mu);

    % Compute true anomaly relative to the computed argument of perigee
    f = theta - argp;

    % Compute radial distance
    r = (h^2 / mu) / (1 + ecc * cos(f));

    % Position and velocity in the **perifocal** frame
    r_vec_eph = [r * cos(f), r * sin(f), 0]; % Position in perifocal frame
    v_vec_eph = [-(mu / h) * sin(f), (mu / h) * (ecc + cos(f)), 0]; % Velocity in perifocal frame

    % Compute rotation matrix for transformation from perifocal to ECI frame
    DCM313 = D313(raan, inc, argp);

    % Transform position and velocity to the **ECI** frame
    r_vec_xyz = DCM313 * r_vec_eph';
    v_vec_xyz = DCM313 * v_vec_eph';

    % Extract Cartesian components
    x = r_vec_xyz(1);
    y = r_vec_xyz(2);
    z = r_vec_xyz(3);
    vx = v_vec_xyz(1);
    vy = v_vec_xyz(2);
    vz = v_vec_xyz(3);
end

% -------------------------------------------------------------------------
% Function: D313
%
% Description:
% Computes the 3-1-3 rotation matrix to transform from the perifocal frame 
% to the ECI frame using RAAN, inclination, and argument of perigee.
%
% Inputs:
%   raan  - Right Ascension of the Ascending Node (RAAN) [radians]
%   inc   - Inclination [radians]
%   argp  - Argument of perigee [radians]
%
% Outputs:
%   D313  - 3x3 Direction Cosine Matrix (DCM)
% -------------------------------------------------------------------------

function [D313] = D313(raan, inc, argp)
    % Construct the 3-1-3 rotation matrix
    D313 = [cos(raan) * cos(argp) - sin(raan) * cos(inc) * sin(argp), ...
           -cos(raan) * sin(argp) - sin(raan) * cos(inc) * cos(argp), sin(raan) * sin(inc); ...
            sin(raan) * cos(argp) + cos(raan) * cos(inc) * sin(argp), ...
           -sin(raan) * sin(argp) + cos(raan) * cos(inc) * cos(argp), -cos(raan) * sin(inc); ...
            sin(inc) * sin(argp), sin(inc) * cos(argp), cos(inc)];
end

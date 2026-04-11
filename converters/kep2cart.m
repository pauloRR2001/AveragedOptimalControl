function [x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,f,mu)
    p = a*(1-ecc^2);
    hmag = sqrt(p*mu);
    rmag = (hmag^2/mu)/(1+ecc*cos(f));
    rVec = [rmag*cos(f), rmag*sin(f), 0];
    vVec = [-(mu/hmag)*sin(f), (mu/hmag)*(ecc+cos(f)), 0];
    DCM313 = D313(raan,inc,argp);
    r = DCM313 *rVec';
    v = DCM313 * vVec';
    x = r(1);
    y = r(2);
    z = r(3);
    vx = v(1);
    vy = v(2);
    vz = v(3);
end

function [D313] = D313(RAAN,inc,argp)
D313 = [cos(RAAN)*cos(argp)-sin(RAAN)*cos(inc)*sin(argp), -cos(RAAN)*sin(argp)-sin(RAAN)*cos(inc)*cos(argp), sin(RAAN)*sin(inc);...
        sin(RAAN)*cos(argp)+cos(RAAN)*cos(inc)*sin(argp), -sin(RAAN)*sin(argp)+cos(RAAN)*cos(inc)*cos(argp), -cos(RAAN)*sin(inc);...
        sin(inc)*sin(argp)                           , sin(inc)*cos(argp)                            , cos(inc)];
end
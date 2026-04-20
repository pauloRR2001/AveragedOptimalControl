function s_cart = KEP2CART(s_kep,mu)
    a = s_kep(1);
    ecc = s_kep(2);
    inc = s_kep(3);
    raan = s_kep(4);
    argp = s_kep(5);
    M = s_kep(6);
    
    E = M2E(M, ecc);
    f = E2f(E,ecc);

    p = a*(1-ecc^2);
    hmag = sqrt(p*mu);
    rmag = (hmag^2/mu)/(1+ecc*cos(f));
    rVec = [rmag*cos(f), rmag*sin(f), 0];
    vVec = [-(mu/hmag)*sin(f), (mu/hmag)*(ecc+cos(f)), 0];

    DCM313 = D313(raan,inc,argp);

    r = DCM313 *rVec';
    v = DCM313 * vVec';

    s_cart = [r; v];
end

function [D313] = D313(RAAN,inc,argp)
D313 = [cos(RAAN)*cos(argp)-sin(RAAN)*cos(inc)*sin(argp), -cos(RAAN)*sin(argp)-sin(RAAN)*cos(inc)*cos(argp), sin(RAAN)*sin(inc);...
        sin(RAAN)*cos(argp)+cos(RAAN)*cos(inc)*sin(argp), -sin(RAAN)*sin(argp)+cos(RAAN)*cos(inc)*cos(argp), -cos(RAAN)*sin(inc);...
        sin(inc)*sin(argp)                           , sin(inc)*cos(argp)                            , cos(inc)];
end
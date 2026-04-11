function Xkep_struct = keplerianPropUnperturbed(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert)
    E0 = 2*atan(((1-ecc0)/(1+ecc0))^0.5 * tan(nu0/2));
    M0 = E0 - ecc0*sin(E0);
    
    Xkep0 = [a0,ecc0,inc0,raan0,argp0,M0];
    Xkep_struct = ode45(@(t,Xkep) dXkepdt(t,Xkep,simT,mu,pert),simT,Xkep0,opt);
end

function dXkepdt = dXkepdt(t,Xkep,simT,mu,pert)
 
    a = Xkep(1);
    ecc = Xkep(2);
    %inc = Xkep(3);
    %raan = Xkep(4);
    %argp = Xkep(5);
    M = Xkep(6);

    %p = a*(1-ecc^2);
    n = sqrt(mu/a^3);
    %h = sqrt(p*mu);
    %b = a*sqrt(1-ecc^2);

    E = meanToEcc(M, ecc);
    nu = 2*atan(((1+ecc)/(1-ecc))^0.5 * tan(E/2));

    %r = p/(1+ecc*cos(nu));

    f0 = [0,0,0,0,0,n]';

    dXkepdt = f0;
end
function ydot = averaged_dynamics(y, params)

    n = 9;

    x = y(1:n);
    lambda = y(n+1:2*n);

    Hfun = @(z) calculate_averaged_hamiltonian(z(1:n), z(n+1:2*n), params);

    z = [x; lambda];

    gradH = numerical_gradient(Hfun, z);

    dH_dx = gradH(1:n);
    dH_dlambda = gradH(n+1:2*n);

    xdot = dH_dlambda;
    lambdadot = -dH_dx;

    ydot = [xdot; lambdadot];
end
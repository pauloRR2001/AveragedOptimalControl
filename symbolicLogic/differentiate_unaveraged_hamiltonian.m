function derivs = differentiate_unaveraged_hamiltonian()

    syms p f g h k Lstate t alph m real
    syms lam_p lam_f lam_g lam_h lam_k lam_L lam_t lam_alph lam_m real
    syms Lq real
    syms mu_g c T_min T_max eps_S RE J2 real

    x = [p; f; g; h; k; Lstate; t; alph; m];
    lambda = [lam_p; lam_f; lam_g; lam_h; lam_k; lam_L; lam_t; lam_alph; lam_m];

    params.mu = mu_g;
    params.c = c;
    params.T_min = T_min;
    params.T_max = T_max;
    params.eps_S = eps_S;
    params.RE = RE;
    params.J2 = J2;

    H = formulate_unaveraged_hamiltonian_symbolic(x, lambda, Lq, params);

    dH_dx = jacobian(H, x).';
    dH_dlambda = jacobian(H, lambda).';

    xdot = dH_dlambda;
    lambdadot = -dH_dx;

    derivs.H = H;
    derivs.x = x;
    derivs.lambda = lambda;
    derivs.Lq = Lq;
    derivs.dH_dx = dH_dx;
    derivs.dH_dlambda = dH_dlambda;
    derivs.xdot = xdot;
    derivs.lambdadot = lambdadot;

    derivs.mu = mu_g;
    derivs.c = c;
    derivs.T_min = T_min;
    derivs.T_max = T_max;
    derivs.eps_S = eps_S;
    derivs.RE = RE;
    derivs.J2 = J2;

end
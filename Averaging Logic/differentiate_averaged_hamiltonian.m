function derivs = differentiate_averaged_hamiltonian(q, use_J2)

    if nargin < 2
        use_J2 = true;
    end

    syms p f g h k L t alph m real
    syms lam_p lam_f lam_g lam_h lam_k lam_L lam_t lam_alph lam_m real
    syms mu_g c T_min T_max epsilon_S eps_u RE J2 real

    x = [p; f; g; h; k; L; t; alph; m];
    lambda = [lam_p; lam_f; lam_g; lam_h; lam_k; lam_L; lam_t; lam_alph; lam_m];

    [z, w] = lgwt(q, -1, 1);
    L_nodes = pi * z;

    params.mu = mu_g;
    params.c = c;
    params.T_min = T_min;
    params.T_max = T_max;
    params.epsilon_S = epsilon_S;
    params.eps_u = eps_u;
    params.RE = RE;
    params.J2 = J2;
    params.use_J2 = use_J2;

    H_avg = calculate_averaged_hamiltonian_symbolic(x, lambda, params, L_nodes, w);

    dH_dx = jacobian(H_avg, x).';
    dH_dlambda = jacobian(H_avg, lambda).';

    derivs.H_avg = H_avg;
    derivs.dH_dx = dH_dx;
    derivs.dH_dlambda = dH_dlambda;
    derivs.xdot = dH_dlambda;
    derivs.lambdadot = -dH_dx;
end
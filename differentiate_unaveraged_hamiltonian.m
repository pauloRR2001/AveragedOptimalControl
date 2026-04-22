function derivs = differentiate_unaveraged_hamiltonian()
    % Symbolic variables
    syms p f g h k L t alph m real
    syms lam_p lam_f lam_g lam_h lam_k lam_L lam_t lam_alph lam_m real
    syms a_r a_th a_h real
    syms mu_g c T_min T_max eps_S real
    
    % Pack state, costate, perturbation
    x = [p; f; g; h; k; L; t; alph; m];
    lambda = [lam_p; lam_f; lam_g; lam_h; lam_k; lam_L; lam_t; lam_alph; lam_m];
    a_pert_rth = [a_r; a_th; a_h];
    
    % Parameter struct
    params.mu = mu_g;
    params.c = c;
    params.T_min = T_min;
    params.T_max = T_max;
    params.eps_S = eps_S;
    
    % Build Hamiltonian
    H = formulate_unaveraged_hamiltonian_symbolic(x, lambda, a_pert_rth, params);
    
    % Derivatives
    dH_dx      = jacobian(H, x).';
    dH_dlambda = jacobian(H, lambda).';
    dH_da      = jacobian(H, a_pert_rth).';
    
    % Pontryagin equations
    xdot      = dH_dlambda;
    lambdadot = -dH_dx;
    
    % Return everything in a struct
    derivs.H = H;
    derivs.x = x;
    derivs.lambda = lambda;
    derivs.a_pert_rth = a_pert_rth;
    derivs.dH_dx = dH_dx;
    derivs.dH_dlambda = dH_dlambda;
    derivs.dH_da = dH_da;
    derivs.xdot = xdot;
    derivs.lambdadot = lambdadot;
end
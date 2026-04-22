function sH = evaluate_sH(x, lambda, L, params)

    x_i = x;
    x_i(6) = L;

    alpha = x_i(8);
    m = x_i(9);

    A = get_MEE_A_matrix(x_i, params.mu);
    B = get_MEE_B_matrix(x_i, params.mu);

    if params.use_J2
        gamma_RTN = J2pert(x_i, params.mu, params.RE, params.J2);
    else
        gamma_RTN = [0; 0; 0];
    end
    
    lam_MEE = lambda(1:6);
    lam_m = lambda(9);
    
    B_trans_lam = B.'*lam_MEE;
    norm_B_trans_lam = norm(B_trans_lam);
    
    if norm_B_trans_lam < 1e-12
        u_hat = [0; 0; 0];
    else
        u_hat = -B_trans_lam/norm_B_trans_lam;
    end

    B_trans_lam = B.'*lam_MEE;
    norm_B_trans_lam = norm(B_trans_lam);

    S = -(params.c/m)*norm_B_trans_lam - lam_m + 1;

    eps_S = params.epsilon_S;

    sigma = 0.5*(1 - S/sqrt(S^2 + eps_S^2));

    T = params.T_min + (params.T_max - params.T_min)*sigma;

    time_costate_term = lambda(7)*alpha;

    cost_term = alpha*(params.T_min/params.c) + ...
                alpha*((params.T_max - params.T_min)/params.c)* ...
                (sigma - eps_S*sqrt(sigma - sigma^2));

    dyn_term = alpha*dot(lam_MEE,A) + ...
               alpha*dot(lam_MEE,B*(u_hat*(T/m) + gamma_RTN));

    mass_costate_term = -alpha*lam_m*(T/params.c);

    H = cost_term + dyn_term + mass_costate_term + time_costate_term;

    p = x_i(1);
    f = x_i(2);
    g = x_i(3);
    
    q = 1 + f*cos(L) + g*sin(L);
    
    if p <= 0 || m <= 0 || q <= 0 || 1 - f^2 - g^2 <= 0
        sH = 1e30;
        return
    end

    r = p/q;

    dt_dL = r^2/sqrt(params.mu*p);

    a = p/(1 - f^2 - g^2);
    P = 2*pi*sqrt(a^3/params.mu);

    s = (2*pi/P)*dt_dL;

    sH = s*H;

end
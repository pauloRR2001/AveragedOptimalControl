function Hsym = formulate_unaveraged_hamiltonian_symbolic(x, lambda, Lq, params)

    p = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    alpha = x(8);
    m = x(9);

    L = Lq;

    lam_MEE = lambda(1:6);
    lam_t   = lambda(7);
    lam_m   = lambda(9);

    q = 1 + f*cos(L) + g*sin(L);
    s = sqrt(1 + h^2 + k^2);

    A = [0;
         0;
         0;
         0;
         0;
         sqrt(params.mu * p) * (q/p)^2];

    B = sqrt(p/params.mu) * ...
        [0,        2*p/q,                                      0; ...
         sin(L),  ((q+1)*cos(L)+f)/q,                         -g*(h*sin(L)-k*cos(L))/q; ...
        -cos(L),  ((q+1)*sin(L)+g)/q,                          f*(h*sin(L)-k*cos(L))/q; ...
         0,        0,                                          (s^2/(2*q))*cos(L); ...
         0,        0,                                          (s^2/(2*q))*sin(L); ...
         0,        0,                                          (h*sin(L)-k*cos(L))/q];

    a_pert_rth = J2pert([p; f; g; h; k; L], params.mu, params.RE, params.J2);

    B_trans_lam = B.' * lam_MEE;
    norm_B_trans_lam = sqrt(B_trans_lam.' * B_trans_lam);

    u_hat = -B_trans_lam / norm_B_trans_lam;

    S = -(params.c / m) * norm_B_trans_lam - lam_m + 1;
    sigma = 0.5 * (1 - S / sqrt(S^2 + params.eps_S^2));

    T = params.T_min + (params.T_max - params.T_min) * sigma;

    cost_rate = alpha * (params.T_min / params.c) + ...
                alpha * ((params.T_max - params.T_min) / params.c) * ...
                (sigma - params.eps_S * sqrt(sigma - sigma^2));

    dyn_term = alpha * (lam_MEE.' * A) + ...
               alpha * (lam_MEE.' * (B * (u_hat * (T / m) + a_pert_rth))) + ...
               alpha * lam_t;

    mass_term = -alpha * lam_m * (T / params.c);

    Hsym = cost_rate + dyn_term + mass_term;

end
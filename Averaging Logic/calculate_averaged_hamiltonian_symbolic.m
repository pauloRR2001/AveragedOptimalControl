function H_avg = calculate_averaged_hamiltonian_symbolic(x, lambda, params, L_nodes, w)

    alph = x(8);
    m = x(9);

    H_avg = sym(0);

    for i = 1:length(L_nodes)
        L_i = L_nodes(i);

        x_i = x;
        x_i(6) = L_i;

        A = get_MEE_A_matrix(x_i, params.mu);
        B = get_MEE_B_matrix(x_i, params.mu);

        if params.use_J2
            gamma_RTN = J2pert(x_i, params.mu, params.RE, params.J2);
        else
            gamma_RTN = sym([0; 0; 0]);
        end

        lam_MEE = lambda(1:6);
        lam_m   = lambda(9);
        lam_t   = lambda(7);

        B_trans_lam = B.' * lam_MEE;
        norm_B_trans_lam = sqrt(B_trans_lam.' * B_trans_lam + params.eps_u^2);
        u_hat = -B_trans_lam / norm_B_trans_lam;

        eps_S = params.epsilon_S;
        S = -(params.c / m) * norm_B_trans_lam - lam_m + 1;
        sigma = 0.5 * (1 - S / sqrt(S^2 + eps_S^2));

        T = params.T_min + (params.T_max - params.T_min) * sigma;

        time_costate_term = lam_t * alph;

        cost_term = alph * (params.T_min / params.c) + ...
                    alph * ((params.T_max - params.T_min) / params.c) * ...
                    (sigma - eps_S * sqrt(sigma - sigma^2));

        dyn_term = alph * (lam_MEE.' * A) + ...
                   alph * (lam_MEE.' * (B * (u_hat * (T / m) + gamma_RTN)));

        mass_costate_term = -alph * lam_m * (T / params.c);

        H_i = cost_term + dyn_term + mass_costate_term + time_costate_term;

        p = x_i(1);
        f = x_i(2);
        g = x_i(3);

        q = 1 + f*cos(L_i) + g*sin(L_i);
        r = p / q;

        dt_dL = r^2 / sqrt(params.mu * p);

        a = p / (1 - f^2 - g^2);
        P = 2*pi*sqrt(a^3 / params.mu);

        s_i = (2*pi / P) * dt_dL;

        H_avg = H_avg + sym(0.5) * w(i) * s_i * H_i;
    end
end
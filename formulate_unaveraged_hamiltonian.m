function H = formulate_unaveraged_hamiltonian(x, lambda, a_pert_rth, params)
    % Extracts the 9 states and 9 costates
    % x = [p, f, g, h, k, L, t, alpha, m]
    % lambda = [lam_p, lam_f, lam_g, lam_h, lam_k, lam_L, lam_t, lam_alpha, lam_m]

    p = x(1); f = x(2); g = x(3); h = x(4); k = x(5); L = x(6);
    alpha = x(8); 
    m = x(9);

    lam_MEE = lambda(1:6);
    lam_t   = lambda(7);
    lam_m   = lambda(9);

    q = 1 + f*cos(L) + g*sin(L);
    s = sqrt(1 + h^2 + k^2);
    
    A = [0; 0; 0; 0; 0; sqrt(params.mu * p) * (q/p)^2];
    
    B = sqrt(p/params.mu) * ...
        [0, 2*p/q, 0; ...
        sin(L), ((q+1)*cos(L)+f)/q, -g*(h*sin(L)-k*cos(L))/q; ...
        -cos(L), ((q+1)*sin(L)+g)/q, f*(h*sin(L)-k*cos(L))/q; ...
        0, 0, (s^2/(2*q))*cos(L); ...
        0, 0, (s^2/(2*q))*sin(L); ...
        0, 0, (h*sin(L)-k*cos(L))/q];

    % Calculate Primer Vector
    B_trans_lam = B' * lam_MEE;
    norm_B_trans_lam = norm(B_trans_lam);
    
    if norm_B_trans_lam > 1e-12
        u_hat = -B_trans_lam / norm_B_trans_lam;
    else
        u_hat = [0; 0; 0]; % Singularity protection
    end

    % Calculate Switching Function (S) and Throttle (sigma)
    S = -(params.c / m) * norm_B_trans_lam - lam_m + 1;
    sigma = 0.5 * (1 - S / sqrt(S^2 + params.eps_S^2));
    
    % Eclipsing function (assuming no shadow for now)
    k_e = 1.0; 
    
    % Actual Thrust magnitude
    T = params.T_min + (params.T_max - params.T_min) * k_e * sigma;

    % Make H
    cost_rate = alpha * (params.T_min / params.c) + ...
                alpha * ((params.T_max - params.T_min) / params.c) * k_e * (sigma - params.eps_S * sqrt(sigma - sigma^2));

    % Trying to make dyn term, needs to be looked over pls
    dyn_term = alpha * dot(lam_MEE, A) + ...
               alpha * dot(lam_MEE, B * (u_hat * (T / m) + a_pert_rth)) + ...
               alpha * lam_t; 

    % mass
    mass_term = -alpha * lam_m * (T / params.c);

    % Total Scalar Hamiltonian
    H = cost_rate + dyn_term + mass_term;
end
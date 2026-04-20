function [H_avg] = calculate_averaged_hamiltonian(x, lambda, params)
    % Extracts states and costates
    % x = [p, f, g, h, k, L, t, alpha, m]
    % lambda = [lam_p, lam_f, lam_g, lam_h, lam_k, lam_L, lam_t, lam_alpha, lam_m]
    
    alpha = x(8); % Time of flight parameter
    m = x(9);     % Mass
   
    %  Set up Quadrature (Brute Force) 
    num_nodes = params.q; 
    
    % Get standard Gauss-Legendre nodes (z) between [-1, 1] and weights (w)
    [z, w] = lgwt(num_nodes, -1, 1); 
    
    % Map z [-1, 1] to True Longitude L [-pi, pi]
    L_nodes = pi * z; 
    
    H_avg = 0; 
    
    % Evaluate Unaveraged Dynamics at each Node 
    for i = 1:num_nodes
        L_i = L_nodes(i);
        
        A = get_MEE_A_matrix(x, L_i, params.mu);
        B = get_MEE_B_matrix(x, L_i, params.mu);
        

        gamma_RTN = get_J2_perturbation(x, L_i, params);
        
        lam_MEE = lambda(1:5);
        B_trans_lam = B' * lam_MEE;
        norm_B_trans_lam = norm(B_trans_lam);
        
        if norm_B_trans_lam > 1e-12
            u_hat = -B_trans_lam / norm_B_trans_lam; % Eq 13
        else
            u_hat = [0; 0; 0]; % Singularity protection
        end
        
        S = -(params.c / m) * norm_B_trans_lam - lambda(9) + 1; 
        
        eps_S = params.epsilon_S;
        sigma = 0.5 * (1 - S / sqrt(S^2 + eps_S^2));
        
        k_e = 1.0; % Eclipsing parameter, set to 1 bc we are ignoring eclipsing for now
        T = params.T_min + (params.T_max - params.T_min) * k_e * sigma;

        cost_term = alpha * (params.T_min / params.c) + ...
                    alpha * ((params.T_max - params.T_min) / params.c) * k_e * (sigma - eps_S * sqrt(sigma - sigma^2));
                    
        dyn_term = alpha * dot(lam_MEE, A) + ...
                   alpha * dot(lam_MEE, B * (u_hat * (T / m) + gamma_RTN));
                   
        mass_costate_term = -alpha * lambda(9) * (T / params.c);
        
        H_i = cost_term + dyn_term + mass_costate_term;
        
        % Calculate s
        % s = (2*pi / P) * (dt/dL)
        p = x(1); f = x(2); g = x(3);
        r = p / (1 + f * cos(L_i) + g * sin(L_i));
        dt_dL = (r^2) / sqrt(params.mu * p);
        
        a = p / (1 - f^2 - g^2); 
        P = 2 * pi * sqrt(a^3 / params.mu);
        
        s_i = (2 * pi / P) * dt_dL;
        
        H_avg = H_avg + 0.5 * w(i) * (s_i * H_i);
    end
end
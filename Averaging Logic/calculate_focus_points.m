function L_stars = calculate_focus_points(x, lambda, params)
    % Extracts the necessary states
    p = x(1); f = x(2); g = x(3); h = x(4); k = x(5); m = x(9);
    lam_MEE = lambda(1:6);
    lam_m = lambda(9);
    beta = 0.05;
    sigma_star_1 = beta;
    sigma_star_2 = 1 - beta;

    % S Values 
    % How to compute S again???
    
    target_S_vals = [S_star_1, S_star_2];

    % Grid 
    L_grid = linspace(-pi, pi, 150); 
    S_vals = zeros(size(L_grid));

    for i = 1:length(L_grid)
        S_vals(i) = evaluate_physical_S(L_grid(i), p, f, g, h, k, m, lam_MEE, lam_m, params);
    end

    L_stars = [];
    options = optimset('Display', 'off'); % Suppress fzero printouts

    for t_idx = 1:2
        S_target = target_S_vals(t_idx);
        residual = S_vals - S_target; % S(L) - S* = 0

        % Find indices where the residual changes sign
        sign_changes = find(residual(1:end-1) .* residual(2:end) < 0);

        for idx = sign_changes
            % IDRK what I'm doing here as of rn, let me look at this again
            % soon
            % L_bracket = [L_grid(idx), L_grid(idx+1)];
            % 
            % % Pinpoint the root
            % root = fzero(@(L) evaluate_physical_S(L, p, f, g, h, k, m, lam_MEE, lam_m, params) - S_target, L_bracket, options);
            % 
            % L_stars = [L_stars, root];
        end
    end

    % Sort the focus points chronologically along the orbit
    L_stars = sort(L_stars);
    
end

function S = evaluate_physical_S(L, p, f, g, h, k, m, lam_MEE, lam_m, params)
    q = 1 + f*cos(L) + g*sin(L);
    s_val = sqrt(1 + h^2 + k^2);
    
    B = sqrt(p/params.mu) * ...
        [0, 2*p/q, 0; ...
        sin(L), ((q+1)*cos(L)+f)/q, -g*(h*sin(L)-k*cos(L))/q; ...
        -cos(L), ((q+1)*sin(L)+g)/q, f*(h*sin(L)-k*cos(L))/q; ...
        0, 0, (s_val^2/(2*q))*cos(L); ...
        0, 0, (s_val^2/(2*q))*sin(L); ...
        0, 0, (h*sin(L)-k*cos(L))/q];
        
    B_trans_lam = B' * lam_MEE;
    S = -(params.c / m) * norm(B_trans_lam) - lam_m + 1;
end
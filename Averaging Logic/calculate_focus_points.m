function L_stars = calculate_focus_points(x, lambda, params)

    beta = 0.05;
    eps_S = params.epsilon_S;

    sigma_star = [beta, 1 - beta];

    target_S_vals = eps_S .* (1 - 2*sigma_star) ./ ...
                    sqrt(1 - (1 - 2*sigma_star).^2);

    L_grid = linspace(-pi, pi, 500);
    S_vals = zeros(size(L_grid));

    for i = 1:length(L_grid)
        S_vals(i) = evaluate_physical_S(L_grid(i), x, lambda, params);
    end

    L_stars = [];

    for j = 1:length(target_S_vals)

        S_target = target_S_vals(j);
        residual = S_vals - S_target;

        sign_changes = find(residual(1:end-1).*residual(2:end) < 0);

        for idx = sign_changes

            L_bracket = [L_grid(idx), L_grid(idx+1)];

            root = fzero(@(L) evaluate_physical_S(L, x, lambda, params) - S_target, ...
                         L_bracket);

            L_stars = [L_stars; root];

        end
    end

    L_stars = sort(L_stars);

    if ~isempty(L_stars)
        keep = [true; abs(diff(L_stars)) > 1e-10];
        L_stars = L_stars(keep);
    end

end

function S = evaluate_physical_S(L, x, lambda, params)

    x_i = x;
    x_i(6) = L;

    m = x_i(9);

    B = get_MEE_B_matrix(x_i, params.mu);

    lam_MEE = lambda(1:6);
    lam_m = lambda(9);

    B_trans_lam = B.' * lam_MEE;

    S = -(params.c/m)*norm(B_trans_lam) - lam_m + 1;

end
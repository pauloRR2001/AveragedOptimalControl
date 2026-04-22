function H_avg = calculate_focused_averaged_hamiltonian(x, lambda, params)

    L_stars = calculate_focus_points(x, lambda, params);
    L_bounds = [-pi; L_stars(:); pi];

    H_avg = 0;

    for j = 1:length(L_bounds)-1

        L_a = L_bounds(j);
        L_b = L_bounds(j+1);

        num_nodes = params.q;
        [z, w] = lgwt(num_nodes, -1, 1);

        for i = 1:num_nodes

            L_i = 0.5*(L_b - L_a)*z(i) + 0.5*(L_b + L_a);
            sH_i = evaluate_sH(x, lambda, L_i, params);

            H_avg = H_avg + (1/(2*pi))*0.5*(L_b - L_a)*w(i)*sH_i;
        end
    end
end
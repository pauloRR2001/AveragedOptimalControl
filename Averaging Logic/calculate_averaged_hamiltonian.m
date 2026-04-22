function H_avg = calculate_averaged_hamiltonian(x, lambda, params)

    num_nodes = params.q;
    [z, w] = lgwt(num_nodes, -1, 1);

    L_nodes = pi*z;

    H_avg = 0;

    for i = 1:num_nodes
        L_i = L_nodes(i);
        sH_i = evaluate_sH(x, lambda, L_i, params);
        H_avg = H_avg + 0.5*w(i)*sH_i;
    end
end
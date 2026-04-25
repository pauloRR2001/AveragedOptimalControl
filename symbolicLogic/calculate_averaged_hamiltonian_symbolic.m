function H_avg = calculate_averaged_hamiltonian_symbolic(x, lambda, params)

[z,w] = lgwt(params.q,-1,1);

H_avg = 0;

for i = 1:length(z)

    Lq = pi*z(i);

    sH_i = sH_integrand_sym_fun( ...
        x, lambda, Lq, ...
        params.mu, params.c, params.T_min, params.T_max, ...
        params.eps_S, params.RE, params.J2);

    H_avg = H_avg + 0.5*w(i)*sH_i;

end

end
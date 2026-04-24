function ydot = averaged_dynamics_symbolic(y, params)

    x = y(1:9);
    lambda = y(10:18);

    [z,w] = lgwt(params.q,-1,1);

    xdot_avg = zeros(9,1);
    dHavg_dx = zeros(9,1);

    for i = 1:length(z)

        Lq = pi*z(i);

        dsHdlam_i = dsHdlam_integrand_sym_fun( ...
            x, lambda, Lq, ...
            params.mu, params.c, params.T_min, params.T_max, ...
            params.eps_S, params.RE, params.J2);

        dsHdx_i = dsHdx_integrand_sym_fun( ...
            x, lambda, Lq, ...
            params.mu, params.c, params.T_min, params.T_max, ...
            params.eps_S, params.RE, params.J2);

        xdot_avg = xdot_avg + 0.5*w(i)*dsHdlam_i;
        dHavg_dx = dHavg_dx + 0.5*w(i)*dsHdx_i;

    end

    lambdadot_avg = -dHavg_dx;

    ydot = [xdot_avg; lambdadot_avg];

end
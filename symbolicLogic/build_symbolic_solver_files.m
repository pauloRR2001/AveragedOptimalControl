clc; clear; close all;

build_symbolic_generated_solve_functions()

function build_symbolic_generated_solve_functions(output_folder)

    if nargin < 1
        output_folder = fullfile(pwd,'generated_symbolic');
    end

    if ~exist(output_folder,'dir')
        mkdir(output_folder);
    end

    derivs = differentiate_unaveraged_hamiltonian();

    x = derivs.x;
    lambda = derivs.lambda;
    Lq = derivs.Lq;

    mu = derivs.mu;
    c = derivs.c;
    T_min = derivs.T_min;
    T_max = derivs.T_max;
    eps_S = derivs.eps_S;
    RE = derivs.RE;
    J2 = derivs.J2;

    p = x(1);
    f = x(2);
    g = x(3);

    q = 1 + f*cos(Lq) + g*sin(Lq);
    r = p/q;
    dt_dL = r^2/sqrt(mu*p);
    a = p/(1 - f^2 - g^2);
    P = 2*pi*sqrt(a^3/mu);
    s = (2*pi/P)*dt_dL;

    sH = s * derivs.H;
    dsHdx = jacobian(sH, x).';
    dsHdlam = jacobian(sH, lambda).';

    d_dsHdlam_dx = jacobian(dsHdlam, x);
    d_dsHdlam_dlam = jacobian(dsHdlam, lambda);
    d_dsHdx_dx = jacobian(dsHdx, x);
    d_dsHdx_dlam = jacobian(dsHdx, lambda);

    vars = {x, lambda, Lq, mu, c, T_min, T_max, eps_S, RE, J2};

    matlabFunction(d_dsHdlam_dx, ...
        'File', fullfile(output_folder,'d_dsHdlam_dx_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', false);

    matlabFunction(d_dsHdlam_dlam, ...
        'File', fullfile(output_folder,'d_dsHdlam_dlam_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', false);

    matlabFunction(d_dsHdx_dx, ...
        'File', fullfile(output_folder,'d_dsHdx_dx_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', false);

    matlabFunction(d_dsHdx_dlam, ...
        'File', fullfile(output_folder,'d_dsHdx_dlam_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', false);

    fprintf('Generated symbolic solver functions in:\n%s\n', output_folder);
    fprintf('  d_dsHdlam_dx_integrand_sym_fun.m\n');
    fprintf('  d_dsHdlam_dlam_integrand_sym_fun.m\n');
    fprintf('  d_dsHdx_dx_integrand_sym_fun.m\n');
    fprintf('  d_dsHdx_dlam_integrand_sym_fun.m\n');

end
clc; clear; close all;

build_symbolic_generated_functions()

function build_symbolic_generated_functions(output_folder)

    if nargin < 1
        output_folder = fullfile(pwd,'generated_symbolic');
    end

    if ~exist(output_folder,'dir')
        mkdir(output_folder);
    end

    delete(fullfile(output_folder,'*.m'))
    clear functions

    derivs = differentiate_unaveraged_hamiltonian();

    x = derivs.x;
    lambda = derivs.lambda;

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

    q = 1 + f*cos(derivs.Lq) + g*sin(derivs.Lq);
    r = p/q;
    dt_dL = r^2/sqrt(mu*p);
    a = p/(1 - f^2 - g^2);
    P = 2*pi*sqrt(a^3/mu);
    s = (2*pi/P)*dt_dL;

    sH = s * derivs.H;
    dsHdx = jacobian(sH, x).';
    dsHdlam = jacobian(sH, lambda).';

    vars = {x, lambda, derivs.Lq, mu, c, T_min, T_max, eps_S, RE, J2};

    matlabFunction(derivs.H, ...
        'File', fullfile(output_folder,'H_unavg_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', true);

    matlabFunction(derivs.xdot, ...
        'File', fullfile(output_folder,'xdot_unavg_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', true);

    matlabFunction(sH, ...
        'File', fullfile(output_folder,'sH_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', true);

    matlabFunction(dsHdx, ...
        'File', fullfile(output_folder,'dsHdx_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', true);

    matlabFunction(dsHdlam, ...
        'File', fullfile(output_folder,'dsHdlam_integrand_sym_fun'), ...
        'Vars', vars, ...
        'Optimize', true);

    fprintf('Generated symbolic functions in:\n%s\n', output_folder);
    fprintf('  H_unavg_sym_fun.m\n');
    fprintf('  xdot_unavg_sym_fun.m\n');
    fprintf('  sH_integrand_sym_fun.m\n');
    fprintf('  dsHdx_integrand_sym_fun.m\n');
    fprintf('  dsHdlam_integrand_sym_fun.m\n');

end
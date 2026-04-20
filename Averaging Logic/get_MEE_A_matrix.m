function A = get_MEE_A_matrix(s, mu)

    p = s(1);
    f = s(2);
    g = s(3);
    L = s(6);

    q = 1 + f*cos(L) + g*sin(L);

    A = [
        0
        0
        0
        0
        0
        sqrt(mu*p)*(q/p)^2
    ];

end
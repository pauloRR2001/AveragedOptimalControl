function Bx = get_MEE_B_matrix(s, mu)

    p = s(1);
    f = s(2);
    g = s(3);
    h = s(4);
    k = s(5);
    L = s(6);

    q = 1 + f*cos(L) + g*sin(L);
    eta2 = 1 + h^2 + k^2;

    Bx = sqrt(p/mu)*[
        0,        2*p/q,                                      0
        sin(L),  ((q+1)*cos(L)+f)/q,                         -g*(h*sin(L)-k*cos(L))/q
       -cos(L),  ((q+1)*sin(L)+g)/q,                          f*(h*sin(L)-k*cos(L))/q
        0,        0,                                          (eta2/(2*q))*cos(L)
        0,        0,                                          (eta2/(2*q))*sin(L)
        0,        0,                                          (h*sin(L)-k*cos(L))/q
    ];

end
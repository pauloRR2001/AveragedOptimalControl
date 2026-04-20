% function a_J2_RTN = J2pert(s_mee,mu)
% 
%     s_kep = MEE2KEP(s_mee);
%     s_cart = KEP2CART(s_kep,mu);
% 
%     r_vec = s_cart(1:3);
%     v_vec = s_cart(4:6);
% 
%     x = r_vec(1);
%     y = r_vec(2);
%     z = r_vec(3);
% 
%     R0 = 6378.1;
%     J2 = 1.0826e-3;
% 
%     r = norm(r_vec);
% 
%     a_J2_cart = -(3*mu*J2*R0^2/(2*r^5))*[
%         (1 - 5*(z/r)^2)*x
%         (1 - 5*(z/r)^2)*y
%         (3 - 5*(z/r)^2)*z
%     ];
% 
%     e_R = r_vec/norm(r_vec);
%     h_vec = cross(r_vec,v_vec);
%     e_N = h_vec/norm(h_vec);
%     e_T = cross(e_N,e_R);
% 
%     C_RTN_I = [
%         e_R.'
%         e_T.'
%         e_N.'
%     ];
% 
%     a_J2_RTN = C_RTN_I*a_J2_cart;
% 
% end
function gamma_RTN = J2pert(s_mee,mu,RE,J2)

    p = s_mee(1);
    f = s_mee(2);
    g = s_mee(3);
    h = s_mee(4);
    k = s_mee(5);
    L = s_mee(6);

    q = 1 + f*cos(L) + g*sin(L);
    r = p/q;

    chi = h*sin(L) - k*cos(L);
    denom = (1 + h^2 + k^2)^2;

    gamma_R = -(3*mu*J2*RE^2/(2*r^4)) * ...
              (1 - 12*chi^2/denom);

    gamma_T = -(12*mu*J2*RE^2/r^4) * ...
              (chi*(k*cos(L) + h*sin(L))/denom);

    gamma_N = -(6*mu*J2*RE^2/r^4) * ...
              ((k*sin(L) - h*cos(L))*(1 - h^2 - k^2)/denom);

    gamma_RTN = [gamma_R; gamma_T; gamma_N];

end
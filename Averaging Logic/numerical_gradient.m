function grad = numerical_gradient(fun, z)

    h = 1e-6;
    grad = zeros(size(z));

    for i = 1:length(z)

        zp = z;
        zm = z;

        dz = h*max(1,abs(z(i)));

        zp(i) = zp(i) + dz;
        zm(i) = zm(i) - dz;

        grad(i) = (fun(zp) - fun(zm))/(2*dz);
    end
end
function f = E2f(E,ecc)
    f = atan2(sin(E) .* sqrt(1 - ecc .^ 2), (cos(E) - ecc));
end
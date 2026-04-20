function [x,w] = lgwt(N,a,b)

    beta = 0.5 ./ sqrt(1 - (2*(1:N-1)).^(-2));

    T = diag(beta,1) + diag(beta,-1);

    [V,D] = eig(T);

    x = diag(D);
    [x,idx] = sort(x);

    V = V(:,idx);

    w = 2*(V(1,:)').^2;

    x = ((b-a)*x + (b+a))/2;
    w = (b-a)*w/2;

end
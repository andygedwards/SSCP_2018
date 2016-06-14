

p = ones(6,1);
A = cio_b(p);

[V,D] = eig(A)
[val,idx] = min(abs(diag(D)))
S_eq = V(:,idx)

A*S_eq


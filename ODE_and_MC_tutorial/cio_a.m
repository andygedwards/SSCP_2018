function dS = func(t,S,p)
% Usage:
% IC = [1,0,0]';
% p = ones(6,1);
% cio_a(0,IC,p)
% [t,S] = ode15s(@cio_a,[0,4],IC,[],p);
% plot(t,S)
% legend('C','I','O')
C = 1;
I = 2;
O = 3; 

k_co = p(1);
k_oc = p(2);
k_oi = p(3);
k_io = p(4);
k_ic = p(5);
k_ci = p(6);

A = zeros(3,3);

A(C, O) = k_co;
A(O, C) = k_oc;
A(O, I) = k_oi;
A(I, O) = k_io;
A(I, C) = k_ic;
A(C, I) = k_ci;

A = A';

for i = 1:length(A)
   A(i,i) = -sum(A(:,i));
end

dS = A*S;


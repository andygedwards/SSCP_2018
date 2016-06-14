

dt = 1;
%dt = 0.1;
t = 0:dt:4;
f = @cio_a;
S = zeros(3,length(t));
IC = [1,0,0]';
p = ones(6,1);
%p(1) = 10;
S(:,1) = IC;
for n = 1:(length(t)-1)
   S(:,n+1) = S(:,n) + dt*f(t, S(:,n), p);
end

plot(t, S);


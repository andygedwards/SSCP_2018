


% use a small time step 
dt = 0.01;
t = 0:dt:2;
% allocated space for states
s = 0*t;

% compute the transition matrix:
p = ones(6,1);
A = cio_b(p);



s(1) = 1; % i.e. start in closed state
for i = 1:(length(t)-1)
    s(i+1) = advance(s(i), A, dt);
end
plot(t,s)

pause

%N = 100;
N = 10000;
S = zeros(length(s),N);

for n = 1:N
  S(1,n)  = 1;
  for i = 1:(length(t)-1)
     S(i+1,n) = advance(S(i,n), A, dt);
  end
end

plot(t,sum(S'==1)/N)

plot(t,sum(S'==1)/N,t,sum(S'==2)/N,t,sum(S'==3)/N)

pause

% do the same using Gillespie


i = 1;
s = 1; % i.e. start in closed state
t = 0;
while t(end)<2
    [s(i+1), dt] = gillespie(s(i), A);
    t(i+1) = t(i) + dt;
    i = i + 1;
end

tfine = 0:0.1:2;
sfine = 0*tfine;
for i = length(t):-1:2
    sfine(find(tfine<t(i))) = s(i-1);
end
plot(t,s,'x',tfine,sfine,'-x')

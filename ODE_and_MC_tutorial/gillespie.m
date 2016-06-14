function [state,dt] = gillespie(state, A);

P = A(:,state); 
P(state) = 0;
lambda = sum(P);

random_number = rand;
dt = -log(random_number)/lambda;


CP = cumsum(P);
CP = CP/CP(end);

random_number = rand;
for i = 1:length(CP)

   if random_number < CP(i)
      state = i;
      return
    end

end
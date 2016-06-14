function state = advance(state, A, dt);

P = A(:,state)*dt; 
P(state) = 0;
CP = cumsum(P);

random_number = rand;
for i = 1:length(CP)

   if random_number < CP(i)
      state = i;
      return
    end

end
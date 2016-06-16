[p,p_names] = systemic_init_parameters();
%p_names(25) 
%i = find(strcmp(p_names, 'SLmin'))
%p(i) = 2.5; %adjust SLmin for isometric twitch

m_names = systemic_monitored_names();
pressure_inds = [5,6,7];
flow_inds = [8,9,10];
Elv_ind = 3;

[init,state_names] = systemic_init_states();


[t,s] = ode15s(@systemic_rhs,[0,4000],init,[],p);
    
pressures = zeros(size(s,1),3);
flows = zeros(size(s,1),3);
Elv = zeros(size(s,1),1);
for i = 1:size(s,1)
   m = systemic_monitor(t(i),s(i,:),p);
   pressures(i,:) = m(pressure_inds);
   flows(i,:) = m(flow_inds);
   Elv(i) = m(Elv_ind);
end

figure(1)
plot(t,pressures)
title('Pressures')

figure(2)
plot(t,flows)
title('Flows')

figure(3) 
plot(t,Elv)
title('LV elastance')


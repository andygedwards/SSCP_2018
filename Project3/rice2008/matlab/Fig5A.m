Slvals = linspace(1.8,2.3,6)
[p,p_names] = rice_model_2008_init_parameters();
p_names(25) 
i = find(strcmp(p_names, 'SLmin'))
p(i) = 2.5; %adjust SLmin for isometric twitch

m_names = rice_model_2008_monitored_names();
force_ind = find(strcmp(m_names,'active'))

[init,state_names] = rice_model_2008_init_states();
SL_ind = find(strcmp(state_names,'SL'))

hold on
for sl  = Slvals
    init(SL_ind) = sl;
   
    [t,s] = ode15s(@rice_model_2008_rhs,[0,1000],init,[],p);
    
    force = zeros(size(s,1),1);
    for i = 1:size(s,1)
        m = rice_model_2008_monitor(t(i),s(i,:),p);
        force(i) = m(force_ind);
    end
    plot(t,force)
end


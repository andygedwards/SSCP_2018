Slvals = linspace(1.8,2.3,6)
[p,p_names] = hybrid_init_parameters();
stretch_id = find(strcmp(p_names,'stretch'));

m_names = hybrid_monitored_names();
force_ind = find(strcmp(m_names,'active'))

[init,state_names] = hybrid_init_states();

hold on
for sl  = Slvals
    %init(SL_ind) = sl;
    s = sl/1.9
    p(stretch_id) = s;
    size(init)
    [t,s] = ode15s(@hybrid_rhs,[0,1000],init,[],p);
    
    force = zeros(size(s,1),1);
    for i = 1:size(s,1)
        m = hybrid_monitor(t(i),s(i,:),p);
        force(i) = m(force_ind);
    end
    plot(t,force)
end


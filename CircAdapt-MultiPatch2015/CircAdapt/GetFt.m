function S= GetFt(Type,Variable,Location)
% function S= Ft(Type,Variable,Location)
% Type=string, Variable=string, Location=string
% e.g. Ft('Node','p','SystArt')

global P;

VarNames=P.(Type).Name;
iLoc=[];
if ischar(Location)
    if strcmp(Location,'All')
        iLoc=1:P.(Type).n;
    else
        iLoc=strmatch(Location,VarNames);
    end
else
    for i=1:length(Location)
        iLoc= [iLoc,strmatch(Location{i},VarNames)];
    end
end
S= P.(Type).(Variable)(:,iLoc);
end


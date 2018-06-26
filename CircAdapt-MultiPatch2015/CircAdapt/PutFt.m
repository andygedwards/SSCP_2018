function PutFt(Type,Variable,Location,Value)
% function PutFt(Type,Variable,Location,Value)
% Type=string, Variable=string, Location=string
% e.g. PutFt('Node','p','SystArt',2000)
% inverse action of GetFt

global P;

Names=P.(Type).Name; %Field names of structure P

iLoc=[]; % column indices
[nt,nVal]=size(Value);

if ischar(Location)
    if strcmp(Location,'All')
        iLoc=1:P.(Type).n;
    else
        iLoc=strmatch(Location,Names);
    end
else
    for i=1:length(Location)
        iLoc= [iLoc,strmatch(Location{i},Names)];
    end
end

nLoc= size(iLoc,2);
if ( nLoc>1 && nVal==1 )
Value=repmat(Value,[1,nLoc]);    
end

if ~isfield(P.(Type),Variable)
    P.(Type).(Variable)= zeros(nt,P.(Type).n);
end
P.(Type).(Variable)(end-nt+1:end,iLoc)=Value;

end


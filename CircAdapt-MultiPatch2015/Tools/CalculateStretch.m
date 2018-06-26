function [SPS,SRS]=CalculateStretch(t,E,iREF,AVO,AVC,invert)

if invert
    E=-E;
end

iREF    = iREF(1);
nCol    = size(E,2);
SRS     = zeros(1,nCol);
SPS     = zeros(1,nCol);

tRes=mean(diff(t));
iAVO=round(AVO/tRes)+1;
iAVC=round(AVC/tRes);

for i = 1:nCol
    
    % systolic pre-stretch (SPS)
    sys_sps = (1:iAVO)+iREF(1)-1; % definition of pre-ejection time interval
    for n = sys_sps(1):sys_sps(end)-1
        if E(n,i)<E(n+1,i),
            dE     = E(n+1,i)-E(n,i);
            SPS(i) = SPS(i)+dE;
        else
        end
    end
    
    % systolic rebound stretch (SRS)
    sys_srs = (1:iAVC)+iREF(1)-1; % definition of systolic time interval
    if min(E(iREF(1):iAVO+iREF(1)-1,i)) < 0
        for n = sys_srs(1):sys_srs(end)-1
            if E(n,i)<E(n+1,i),
                dE     = E(n+1,i)-E(n,i);
                SRS(i) = SRS(i)+dE;
            else
            end
        end
    else
    end
end

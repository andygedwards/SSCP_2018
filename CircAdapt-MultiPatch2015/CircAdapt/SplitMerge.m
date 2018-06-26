function SplitMerge(PatchName,nSM)
% function SplitMerge(PatchName,nSM);
% SplitMerge splits or merges patches within a wall
% PatchName is string
% nSM= number of patches after splitting,
% if nSM<0, - number of merging patches, starting from named patch
% if more merging than possible, than nSM merging is maximized.
% Example: PatchName='Lv1', nSM=3-> LV patch split in 3 patches
% followed by: PatchName='Lv1', nSM=-2 -> LV patch 2,3 merged
% All variables are copied, wall volume and area is subdivided/averaged

global P
% PatchName='Lv1'; nSM=-2;

WallName=PatchName(isletter(PatchName)); % wall identity
nP=GetFt('Wall','nPatch',WallName); % current number of patches in wall
FieldPatch= fieldnames(P.Patch); % all field names in P.Patch
nfield=length(FieldPatch);
FieldAdapt= fieldnames(P.Patch.Adapt); % all field names in P.Patch
nAdapt=length(FieldAdapt);
iP=strmatch(PatchName,P.Patch.Name); %location to insert Patch.records
if nSM>1
    PutFt('Wall','nPatch',WallName,nP+nSM-1); % new value of nPatch
    P.Patch.n= P.Patch.n+nSM-1; % new number of patches
    RgP=iP+(0:nSM-1); % indices of patch group
    for i=1:nAdapt % copying adaptation target values
        Var= P.Patch.Adapt.(FieldAdapt{i});
        nVar=size(Var,2);
        if nVar>1
            Var=Var(:,[1:iP-1,repmat(iP,[1,nSM]),iP+1:end]);
            P.Patch.Adapt.(FieldAdapt{i})=Var;
        end
    end
    for i=1:nfield % copy fields
        Var= P.Patch.(FieldPatch{i});
        nVar=size(Var,2);
        if nVar>1
            Var=Var(:,[1:iP-1,repmat(iP,[1,nSM]),iP+1:end]);
            P.Patch.(FieldPatch{i})=Var;
        end
    end
    % sub division in smaller part for volumes and areas
    P.Patch.VWall(RgP)=P.Patch.VWall(RgP)/nSM;
    P.Patch.AmRef(RgP)=P.Patch.AmRef(RgP)/nSM;
    P.Patch.DADT(:,RgP)=P.Patch.DADT(:,RgP)/nSM;
    P.Patch.Am0(:,RgP)=P.Patch.Am0(:,RgP)/nSM;
    P.Patch.Am(:,RgP)=P.Patch.Am(:,RgP)/nSM;
    P.Patch.dT(RgP)=P.Patch.dT(RgP)/nSM;
else % Merge is inverse of Split
    iPMax=GetFt('Wall','iPatch',WallName)+nP-1;
    nMerge=min(iPMax-iP+1,-nSM); %number of patches to merge
    PutFt('Wall','nPatch',WallName,nP-nMerge+1); % new value of nPatch
    P.Patch.n= P.Patch.n-nMerge+1; % new number of patches
    RgP=iP+(0:nMerge-1); % index patches to merge
    for i=1:nfield
        Var= P.Patch.(FieldPatch{i});
        nVar=size(Var,2);
        if nVar>1
            if isnumeric(Var)
                Var(:,iP)= mean(Var(:,RgP),2);
            end
            Var=Var(:,[1:iP,iP+nMerge:end]);
            P.Patch.(FieldPatch{i})=Var;
        end
    end
    for i=1:nAdapt
        Var= P.Patch.Adapt.(FieldAdapt{i});
        nVar=size(Var,2);
        if nVar>1
            if isnumeric(Var)
                Var(:,iP)= mean(Var(:,RgP),2);
            end
            Var=Var(:,[1:iP,iP+nMerge:end]);
            P.Patch.Adapt.(FieldAdapt{i})=Var;
        end
    end
    % volumes and areas have to be summed instead of averaged
    P.Patch.VWall(iP)=P.Patch.VWall(iP)*nMerge;
    P.Patch.AmRef(iP)=P.Patch.AmRef(iP)*nMerge;
    P.Patch.DADT(:,iP)=P.Patch.DADT(:,iP)*nMerge;
    P.Patch.Am0(:,iP)=P.Patch.Am0(:,iP)*nMerge;
    P.Patch.Am(:,iP)=P.Patch.Am(:,iP)*nMerge;
    P.Patch.dT(iP)=P.Patch.dT(iP)*nMerge;
end

n2iPatch; % finish

end

function n2iPatch
global P
P.Wall.iPatch=[1,1+cumsum(P.Wall.nPatch(1:end-1))]; %renumber patches
NamePatch={}; % new patch names
for i=1:P.Wall.n
    n=P.Wall.nPatch(i);
    for j=1:n
        NamePatch=[NamePatch,[P.Wall.Name{i},num2str(j)]];
    end
end
P.Patch.Name=NamePatch;
P2SVar; % new set of state variables
end


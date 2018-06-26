function strTot=MapStructure(Struct,varargin)
%function MapStructure(Struct,varargin);
%Commonly used as function MapStructure(Struct)
%Nested search in structure tree to map sizes of branch arrays within structure
%Single records (end branches) are not printed
%The varargin facility is needed for the specific nesting in this procedure

%INPUT
% Struct=stucture to map
% varargin{1}=nr, current branch number
% varargin{2}=depth, depth in structure tree
% varargin{3}=string array of output
%OUTPUT
% cell array of strings, containing output to print
%Theo Arts
%==== initialization of Root of structure

LenMin=0; %determines whether compact (LenMin=1) or extended map (LenMin=0)
if nargin<4; % initialization of Root of structure
   nr=1; depth=0; strTot={};
   if isfield(Struct,'Name');
      strTot=[strTot;{Struct.Name}];
      %disp(Struct.Name);
   end
else
   nr=varargin{1}; depth=varargin{2}; strTot=varargin{3};
end

if ~isstruct(Struct); return; end; % Escape if not a structure
Fn1=fieldnames(Struct);
n1=length(Fn1);

for i1=1:n1,
  Twig=getfield(Struct,{1},Fn1{i1});
  Len=length(Twig);
  if (~ischar(Twig)) && (Len>LenMin);
%       str=[repmat(' ',[1,4*depth]), num2str(nr),'  ',Fn1{i1},'(',num2str(Len),')'];
      str=[repmat(' ',[1,4*depth]), num2str(nr),'  ',Fn1{i1},'(',num2str(Len),')'];
      strTot=[strTot;{str}];
      %disp(str);
      for i2=1:Len
          strTot=MapStructure(Twig(i2),i2,depth+1,strTot);
      end
  end
end

return


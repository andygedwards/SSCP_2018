function OutNew=SteadyStateP
% function OutNew=SteadyStateP
% OutNew= new estimate vector of steady state parameters, to be used 
% as start condition for next heart beat


global P

% read input-output record, take logarithm for relative change
In  = log(P.Adapt.In );
Out = log(P.Adapt.Out);
[nx,np]=size(In); % [number of samples, number of parameters]
yRef= Out(end,:);

nxn= ceil(0.5*np);
nMin= max(2,nx-nxn); % number of samples used for prediction
Rgn=(nMin:nx)'; %indices used for prediction
nxn=length(Rgn);

if nxn>1
    X=In(Rgn,:); Y=Out(Rgn,:); %U=Out(Rgn-1,:);
    nt=size(X,1); Col1=ones(nt,1);
    YmX =Y-X;
    M1 = [YmX ,Col1];
    N1 = pinv(M1)*(Y-Col1*yRef);
    dy  =  N1(end,:);
    ymx=YmX(end,:); %last instationarity difference
%     ymx=X(end,:)-Y(end-1,:); %last dy
    fac=sqrt( sum(ymx.^2)/sum(dy.^2+1e-10) );
    dy= min(fac,1.2)*dy; % correction amplitude less than last difference
    yNew=yRef+dy;
else
    yNew=yRef;
end
    
OutNew=exp(yNew); % prediction after exponential

end
%     xm=mean(X); ym=mean(Y);
%     X1=X-Col1*xm; Y1=Y-Col1*ym;
%     DXDY= pinv(Y1)*X1;
%     yNew= (yRef-ym)*DXDY + yRef;

%     YRef= repmat(yRef,[nxn,1]);
%     Fac=10; w0=1;
%     while nx1>1 && Fac>1.0
%         Rg=(nxn-nx1+1:nxn)';
%         X1  = X(Rg,:)-YRef(Rg,:);
%         Y1  = Y(Rg,:)-YRef(Rg,:);
%         U1  = U(Rg,:)-YRef(Rg,:);
%         w= w0*sqrt(mean((Y1(:)-X1(:)).^2));
%         Col1= w*ones(nx1,1);
%         
%         Method=1;
%         if Method==1 % no history dependency
%             YmX =Y1-X1;
%             M1  = [YmX ,Col1];
%             N1 = pinv(M1)*Y1;
%             y0  =  w*N1(end,:);
%         end;
%         if Method==2 % history dependency
%             YmXU= [Y1-X1,Y1-U1];
%             M1  = [YmXU,Col1];
%             N1 = pinv(M1)*Y1;
%             y0 =  w*N1(end,:);
%         end;
% 
%         Fac=norm(y0)/(ErrY0*sqrt(np));
%         nx1=nx1-1;
%     end

%     disp(num2str([nxn,nx1+1]));
%     disp(round(1000*[y0;Y(end,:)-X(end,:)]));
% % % % end
% % % % % disp(['Number of used samples for SteadyState: ',num2str(nx1)]);
% % % % % disp(num2str(y0));
% % % % y1= yRef+y0;
% % % % 


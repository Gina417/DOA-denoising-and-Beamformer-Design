function y_hat_LCMV= LCMV(matX,theta_s_hat,theta_i_hat,as,ai,sigma_o,mu)

N=length(matX(:,1)); %number of sensors(N=10)     
L=length(matX(1,:)); %length of time(L=2000) 

P=eye(N)/(sigma_o^2);

y_hat_LCMV=[];

w=as/N;
    
g_con=[1;10^-4];

for m=1:L
  
    C=[as(:,m) ai(:,m)];
    g=(1/mu)*(P*matX(:,m))/(1+(1/mu)*matX(:,m)'*P*matX(:,m));
    P=(1/mu)*P-(1/mu)*(g*matX(:,m)'*P);    % P=inv(R);
    w=P*C*inv(C'*P*C)*g_con;
    y_hat_LCMV(:,m)=w'*matX(:,m);
    
end
end
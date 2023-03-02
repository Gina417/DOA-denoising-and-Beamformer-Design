function y_hat_MVDR = MVDR(matX,theta_s_hat,as,sigma_o,mu)

N=length(matX(:,1)); %number of sensors(N=10)     
L=length(matX(1,:)); %length of time(L=2000)      

P=eye(N)/(sigma_o^2);

y_hat_MVDR=[];

w=as/N;
    
for m=1:L

    g=(1/mu)*(P*matX(:,m))/(1+(1/mu)*matX(:,m)'*P*matX(:,m));
    P=(1/mu)*P-(1/mu)*(g*matX(:,m)'*P);    % P=inv(R);
    w=(P*as(:,m))/(as(:,m)'*P*as(:,m));
    y_hat_MVDR(:,m)=w'*matX(:,m);
end
end





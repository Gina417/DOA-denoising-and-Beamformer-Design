function s_t_hat = my_beamformer(matX,theta_s_hat,theta_i_hat,as,ai_mi,ai,ai_pl,sigma_o,mu)

N=length(matX(:,1)); %number of sensors(N=10)     
L=length(matX(1,:)); %length of time(L=2000) 

P=eye(N)/(sigma_o^2);

s_t_hat=[];

w=as/N;


for m=1:L
  
    C=[as(:,m) ai_mi(:,m) ai(:,m) ai_pl(:,m)];
    g=(1/mu)*(P*matX(:,m))/(1+(1/mu)*matX(:,m)'*P*matX(:,m));
    P=(1/mu)*P-(1/mu)*(g*matX(:,m)'*P); % P=inv(R);

    if (m>=670 && m<=1172) 
        g_con=[1;10^-4;10^-4;10^-4];
    else
        g_con=[1;0;0;0];
    end

    w=P*C*inv(C'*P*C)*g_con;
    s_t_hat(:,m)=w'*matX(:,m);
end
end
function y_hat_as = array_steering(matX,theta_s_hat,as)

N=length(matX(:,1)); %number of sensors(N=10)_
L=length(matX(1,:)); %length of time(L=2000)

y_hat_as=[];

for m=1:L
   
        w(:,m)=as(:,m)/N;
        y_hat_as(:,m)=w(:,m)'*matX(:,m);
    
end

end
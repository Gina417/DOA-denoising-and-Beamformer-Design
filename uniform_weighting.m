function y_hat_uw = uniform_weighting(matX)

N=length(matX(:,1)); %number of sensors(N=10)
L=length(matX(1,:)); %length of time(L=2000)

w=ones(N,1)/N;

for m=1:L

    y_hat_uw(:,m)=w'*matX(:,m);
  
end

end
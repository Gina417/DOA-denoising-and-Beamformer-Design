function y = EMD(x,t,thr)

y=[];
k=1;
w=x;
x0=x;

while 1
    while 1
        wlp=islocalmax(w);
        lpp=t(wlp); %the position of local peaks
        lp=w(wlp); %local peaks

        lpp=[t(1) lpp t(end)];
        lp=[w(1) lp w(end)];

        e_max=spline(lpp,lp,t); %connect local peaks
     
        wld=islocalmin(w);
        ldp=t(wld); %the position of local dips
        ld=w(wld); %local dips
        
        ldp=[t(1) ldp t(end)];
        ld=[w(1) ld w(end)];

        e_min=spline(ldp,ld,t); %connect local dips

        m=(1/2)*(e_max+e_min);
        d=w-m;

        wlp_d=islocalmax(d);
        lpp_d=t(wlp_d); %the position of local peaks of hk
        lp_d=d(wlp_d);

        lpp_d=[t(1) lpp_d t(end)];
        lp_d=[d(1) lp_d d(end)];

        wld_d=islocalmin(d);
        ldp_d=t(wld_d); %the position of local dips of hk
        ld_d=d(wld_d);

        d_e_max=spline(lpp_d,lp_d,t);
        d_e_min=spline(ldp_d,ld_d,t);

        u=(1/2)*(d_e_max+d_e_min);
        k=k+1;
        
        if (all(lp_d>0) && all(ld_d<0) && all(abs(u)<thr) || k>1000)
           y=[y;d];
           x0 = x0 - d;
           break
        else
            w=d;
        end
        
    end
      
    wlp_x0=islocalmax(x0);
    lp_x0=x0(wlp_x0);

    wld_x0=islocalmin(x0);
    ld_x0=x0(wld_x0);

   num_lp_x0=length(lp_x0);
   num_ld_x0=length(ld_x0);

   if (num_lp_x0+num_ld_x0 <= 3)
       break
   else
       w=x0;
   end
end

y = [y;x0];
end
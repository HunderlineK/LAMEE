%Empirical Solution
m=1
n=10000
for p=m:n
   %empirical
   h_fg(p) = 2.501 - 0.0024 * (Tsol_in_(p));
   cp_sol(p) = HeatCapacity(salt_{p}, (Tsol_in_(p)), Xsol_in_(p));
   NTU_m(p) = NTU_(p) * NTU_ratio_(p);
   
%    eff_r(p) = eff_m_emp(p)/eff_s_emp(p);
%    if eff_r(p)<= 0.3 && eff_r(p)>=0
%        eff_r(p) = 0.3
%    elseif eff_r(p)>= -0.1 && eff_r(p)<=0
%        eff_r(p) = -0.1
%    end

   cp_sol(p) = HeatCapacity(salt_{p}, (Tsol_in_(p)), Xsol_in_(p));
   NTU_m(p) = NTU_(p) * NTU_ratio_(p);
   
   eff_r(p) = (1-exp(-NTU_m(p)))/(1-exp(-NTU_(p)));
   Cr_e_emp(p) = Cr_(p) * (1 + h_star(p) * eff_r(p));
   eff_s_emp(p) = eff(NTU_(p), Cr_e_emp(p));
   dTair_emp(p) = eff_s_emp(p)*(Tsol_in_(p) - Tair_in_(p));
   
   mse_in(p) = 0.056 * Wsol_in_(p) * Cr_(p) * h_fg(p) * (1+1/(h_star(p) * eff_r(p)));
   mse_IWD(p) = 0.056 * (Wair_in_(p)-Wsol_in_(p)) * Cr_(p) * h_fg(p) * (1 + 1 /(h_star(p) * eff_r(p)));   
   f = @(x) (x - m_star_e_de(mse_in(p), mse_IWD(p), NTU_m(p),x));
   
   fval = 1;
   trial = 1;
   options = optimset('Display','off');
   while fval>=0.01 && trial <= 150
       if trial == 1
           [m_star_e_decoupled(p),fval] = lsqnonlin(f, 25,-100,100,options);
       else
           [m_star_e_decoupled(p),fval] = lsqnonlin(f, random('unif',-100,100),-1e3,1e3,options);
       end
       trial = trial+1;
   end
   if fval>0.01
       failed_emp(p)=1;
   else
       failed_emp(p)=0;
   end
   mstar{p}=[p, mse_in(p),mse_IWD(p),m_star_e_decoupled(p)];
   
   eff_m_emp(p) = eff(NTU_m(p), m_star_e_decoupled(p));
   dWair_emp(p) = eff_m_emp(p)*(Wsol_in_(p) - Wair_in_(p));
   dHair_emp(p) = dTair_emp(p) + h_fg(p) * dWair_emp(p);

   Tair_out_emp (p) = Tair_in_(p) + dTair_emp(p);
   Wair_out_emp (p) = Wair_in_(p) + dWair_emp(p);
  
end

err = eff_T_num - eff_T_emp;
   figure
   plot(dHair_emp(m:n),dHair_num(m:n),'.')
   hold on
   grid on
   plot([-100,100],[-100,100],'--r')
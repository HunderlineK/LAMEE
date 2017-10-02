% Sanity test

clc
z=0;
for i=22:22 
   %empirical
   h_fg(i) = min(2.501 - 0.0024 * Tsol_in_(i),0.5);
   cp_sol(i) = HeatCapacity(salt_, Tsol_in_(i), Xsol_in_(i));
   NTU_m(i) = NTU_(i) * NTU_ratio_(i);
   
   eff_r(i) = (1-exp(-NTU_m(i)))/(1-exp(-NTU_(i)));
   h_star(i) = h_fg(i) * ( Wsol_in_(i) - Wair_in_(i))/(Tsol_in_(i) - Tair_in_(i));
   Cr_e_emp(i) = Cr_(i) * (1 + h_star(i) * eff_r(i));
   eff_s_emp(i) = eff(NTU_(i), Cr_e_emp(i));
   dTair_emp(i) = eff_s_emp(i)*(Tsol_in_(i) - Tair_in_(i));
   
   mse_in(i) = 0.058 * Wsol_in_(i) * Cr_(i) * h_fg(i) * (1+1/(h_star(i) * eff_r(i)));
   mse_IWD(i) = 0.058 * (Wair_in_(i)-Wsol_in_(i)) * Cr_(i) * h_fg(i) * (1 + 1 /(h_star(i) * eff_r(i)));
   f = @(x) (x - m_star_e_de(mse_in(i), mse_IWD(i), NTU_m(i),x));
   options = optimset('Display','off');
   m_star_e_decoupled(i) = fsolve(f, 0.5,options);
   eff_m_emp(i) = eff(NTU_m(i), m_star_e_decoupled(i));
   dWair_emp(i) = eff_m_emp(i)*(Wsol_in_(i) - Wair_in_(i));
   

   dHair_emp(i) = dTair_emp(i) + h_fg(i) * dWair_emp(i);
   
   % Numerical

   x01 =  40; %+ dTair_emp(i); %Tair_in_(i) + dTair_emp(i);
   x02 =  30; %+ dWair_emp(i); %Wair_in_(i) + dWair_emp(i);
   options = optimset('Display','iter','MaxFunEvals',1000);
   [x,fval,eflag] = lsqnonlin(@(x) LAMEE_ode(salt_, NTU_(p) , NTU_ratio_(p), ...
   Cr_(p),...
   x(1),x(2), Tsol_in_(p), Xsol_in_(p))... 
   -[ Tair_in_(p), Wair_in_(p)],[x01 ,  x02],[2 2], [60 60],options);

   if abs(eflag*eflag')>1
       p
       abs(eflag*eflag')
   end
   
   Tsol_out(i) = gTsol(length(gTsol));
   Tair_out(i) = real(gTair(1));
   Cr_e_num(i) = -(Tsol_out(i) - Tsol_in_(i))/(Tair_out(i) - Tair_in_(i));
   dTair_num(i) = Tair_out(i) - Tair_in_(i);
   %eff_s_num(i) = (Tair_out(i) - Tair_in(i))/(Tsol_in(i) - Tair_in(i));
   
   Wair_out(i) = real(gWair(1));
   dWair_num(i) = Wair_out(i) - Wair_in_(i);
   %eff_m_num(i) = (Wair_out(i) - Wair_in(i))/(Wsol_in(i) - Wair_in(i));
   
   %eff_T_num(i) = (eff_s_num(i) + h_star(i) * eff_m_num(i))/(1 + h_star(i));
   dHair_num(i) = dTair_num(i) + h_fg(i) * dWair_num(i);
   
   Table{i,1} = gTau;
   Table{i,2} = gTair;
   Table{i,3} = gWair;
   Table{i,4} = gTsol;
   Table{i,5} = gXsol;
%end
end
z
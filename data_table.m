clc
warning('off');
options = optimset('Display','off','FinDiffType','central','MaxFunEvals',1000, 'MaxIter',1000);

global gTau gTair gWair gTsol gXsol;

clear dTair_emp
clear dTair_num
clear dWair_emp
clear dWair_num

Table = cell(91,5);

for i=1:1000
 
%    NTU_i = NTU(i);
%    NTU_ratio_i = NTU_ratio(i);
%    salt_i = Salt(i);
%    Tsol_in_i = Tsol_in(i);
%    Xsol_in_i = Xsol_in(i);
%    Wsol_in_i = W(VP(salt_i, Tsol_in_i, Xsol_in_i));
%    Wair_in_i = Wair_in(i);
%    Tair_in_i = Tair_in(i);
   
   Cr(i)= random('unif',0.001,0.999);
   NTU_i = random('unif',1,8);
   NTU_ratio_i = random('unif',0.01,0.99);
   salt_i = 'LiCl';
   Tsol_in_i = random('unif',10,35);
   Xsol_in_i = random('unif',0.01,35);
   Wsol_in_i = W(VP(salt_i, Tsol_in_i, Xsol_in_i));
   Tair_in_i = random('unif',10,45);
   Wsat_air_in_i = W(VP(salt_i, Tair_in_i, 0.01));
   Wair_in_i = random('unif',0.01, Wsat_air_in_i);
   
   
   %empirical
   h_fg = 2.501 - 0.0024 * Tsol_in_i;
   cp_sol = HeatCapacity(salt_i, Tsol_in_i, Xsol_in_i);
   NTU_m_i = NTU_i * NTU_ratio_i;
   
   eff_r = (1-exp(-NTU_m_i))/(1-exp(-NTU_i));
   h_star(i) = h_fg * ( Wsol_in_i - Wair_in_i)/(Tsol_in_i - Tair_in_i);
   Cr_e_emp(i) = Cr(i) * (1 + h_star(i) * eff_r);
   eff_s_emp(i) = eff(NTU_i, Cr_e_emp(i));
   dTair_emp(i) = eff_s_emp(i)*(Tsol_in_i - Tair_in_i);
   
   mse_in = 0.058 * Wsol_in_i * Cr(i) * h_fg * (1+1/(h_star(i) * eff_r));
   mse_IWD = 0.058 * (Wair_in_i-Wsol_in_i) * Cr(i) * h_fg * (1 + 1 /(h_star(i) * eff_r));
   f = @(x) (x - m_star_e_de(mse_in, mse_IWD, NTU_m_i,x));
   options = optimset('Display','off');
   m_star_e_decoupled(i) = fsolve(f, 0.5,options);
   eff_m_emp(i) = eff(NTU_m_i, m_star_e_decoupled(i));
   dWair_emp(i) = eff_m_emp(i)*(Wsol_in_i - Wair_in_i);
   
   % Numerical
   x = fsolve(@(x) LAMEE_ode_BDF(salt_i, NTU_i , NTU_ratio_i, ...
       Cr(i),...
       x(1),x(2), Tsol_in_i, Xsol_in_i)... 
       -[ Tair_in_i, Wair_in_i],[Tsol_in_i + 0.5 ,  Wsol_in_i+0.5],options)
   
   Tsol_out_i = gTsol(length(gTsol));
   Tair_out_i = real(gTair(1));
   Cr_e_num(i) = -(Tsol_out_i - Tsol_in_i)/(Tair_out_i - Tair_in_i);
   dTair_num(i) = Tair_out_i - Tair_in_i;
   eff_s_num(i) = (Tair_out_i - Tair_in_i)/(Tsol_in_i - Tair_in_i);
   
   Wair_out_i = real(gWair(1));
   dWair_num(i) = Wair_out_i - Wair_in_i;
   eff_m_num(i) = (Wair_out_i - Wair_in_i)/(Wsol_in_i - Wair_in_i);
   
   Table{i,1} = gTau;
   Table{i,2} = gTair;
   Table{i,3} = gWair;
   Table{i,4} = gTsol;
   Table{i,5} = gXsol;

end
    
%    figure
%    plot(dTair_num,'--r')
%    hold on
%    plot(dTair_emp)

%    figure
%    plot(dTair_num,dTair_emp,'o')
   
   figure
   plot(dTair_num,dTair_emp,'.')
   grid on
   
   figure
   plot(eff_s_num,eff_s_emp,'.')
   grid on
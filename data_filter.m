%Filter Data from Saturation

clear Cr_filtered;
clear NTU_filtered;
clear NTU_m_filtered;
clear h_star_filtered;

clear dHair_emp_filtered;
clear dHair_num_filtered;
clear dTair_emp_filtered;
clear dTair_num_filtered;
clear dWair_emp_filtered;
clear dWair_num_filtered;

clear Tair_in_filtered;
clear Wair_in_filtered;
clear Tair_out_filtered;
clear Wair_out_filtered;

clear Tsol_in_filtered;
clear Xsol_in_filtered;
clear Wsol_in_filtered;
clear Tsol_out_filtered;
clear Xsol_out_filtered;
clear Wsol_out_filtered;

clear eff_T_emp_filtered;
clear eff_T_num_filtered;
clear eff_m_emp_filtered;
clear eff_m_num_filtered;
clear eff_s_emp_filtered;
clear eff_s_num_filtered;

q = 1;
for p = m:n
    %1/(h_star(p)+1) <= 0.7 && 1/(h_star(p)+1)>=0.1 && h_star(p)>=1 && failed_emp(p)==0 && failed_num(p)==0
    if  1/(h_star(p)+1) <= 10 && 1/(h_star(p)+1)>=-5
        Cr_filtered(q) = Cr_(p);
        NTU_filtered(q) = NTU_(p);
        NTU_m_filtered(q) =  NTU_(p) * NTU_ratio_(p);
        
        h_star_filtered(q) = h_star(p);
        
    	Tair_in_filtered(q) = Tair_in_(p);
        Wair_in_filtered(q) = Wair_in_(p);
        Hair_in_filtered(q) = Tair_in_filtered(q) + 2.4 * Wair_in_filtered(q);
        Tair_out_filtered(q) = Tair_out(p);
        Wair_out_filtered(q) = Wair_out(p);
        Hair_out_filtered(q) = Tair_out_filtered(q) + 2.4 * Wair_out_filtered(q);
        
        Tsol_in_filtered(q) = Tsol_in_(p);
        Xsol_in_filtered(q) = Xsol_in_(p);
        Wsol_in_filtered(q) = W(VP(salt_{p}, Tsol_in_filtered(q), Xsol_in_filtered(q)));
        Hsol_in_filtered(q) = Tsol_in_filtered(q) + 2.4 * Wsol_in_filtered(q);
        Tsol_out_filtered(q) = Tsol_out(p);
        %Xsol_out_filtered(q) = Table{p,5}(length(Table{p,5}));
        %Wsol_out_filtered(q) = W(VP('LiCl', Tsol_out_filtered(q), Xsol_out_filtered(q)));
        %Hsol_out_filtered(q) = Tsol_out_filtered(q) + 2.4 * Wsol_out_filtered(q);
        
        mstar_filtered(q) = Cr_(p) * HeatCapacity('LiCl', Tsol_in_filtered(q), Xsol_in_filtered(q));
        
        eff_s_Cr_filtered(q) = eff(NTU_(p),Cr_(p));
        eff_m_mstar_filtered(q) = eff(NTU_m(p),mstar_filtered(q));
        eff_T_filtered(q) = (eff_s_Cr_filtered(q) + h_star(p)*eff_m_mstar_filtered(q))/(1+ h_star(p));
        dHair_filtered(q) = eff_T_filtered(q) * (Hsol_in_filtered(q) - Hair_in_filtered(q));
        
        dTair_emp_filtered(q) = dTair_emp(p);
        dTair_num_filtered(q) = dTair_num(p);
        eff_s_emp_filtered(q) = eff_s_emp(p);
        eff_s_num_filtered(q) = eff_s_num(p);
        
        dWair_emp_filtered(q) = dWair_emp(p);
        dWair_num_filtered(q) = dWair_num(p);
        eff_m_emp_filtered(q) = eff_m_emp(p);
        eff_m_num_filtered(q) = eff_m_num(p);
        
        dHair_emp_filtered(q) = dHair_emp(p);
        dHair_num_filtered(q) = dHair_num(p);
        eff_T_emp_filtered(q) = eff_T_emp(p);
        eff_T_num_filtered(q) = eff_T_num(p);
        
        
        eff_r_filtered(q) = eff_r(p);
        
        q = q + 1;
    end
end
q

figure; plot(h_star_filtered, dTair_emp_filtered-dTair_num_filtered,'.');
%figure; plot(h_star_filtered,(eff_s_emp_filtered-eff_s_num_filtered),'.'); grid;

figure; plot(h_star_filtered, dWair_emp_filtered-dWair_num_filtered,'.');
%figure; plot(h_star_filtered,(eff_m_emp_filtered-eff_m_num_filtered),'.'); grid;


figure; plot(h_star_filtered,(dHair_emp_filtered-dHair_num_filtered).*Cr_filtered,'.'); grid;
figure; plot(h_star_filtered,(eff_T_emp_filtered-eff_T_num_filtered),'.'); grid;

figure; plot(Tair_in_filtered,Wair_in_filtered,'.');
hold; plot(Tsol_in_filtered,Wsol_in_filtered,'xr');

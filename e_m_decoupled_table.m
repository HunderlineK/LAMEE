clc

for p = 0;%[0.5 1 1.5 2.5 5 10]
    mse_IWD =  p; 
    hold on;
    is_dashed = false;
    for q=[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5 3 4 5]
        mse_in = q
        NTU(1) = 0;
        m_star_e_decoupled(1) = mse_in;
        e(1) = 0;
        options = optimset('Display','off','Algorithm','trust-region-reflective');
        for k = 2:51
            NTU(k) = (k-1)*0.1;
            f = @(x) (x - m_star_e_de(mse_in,mse_IWD,NTU(k),x));
            m_star_e_decoupled(k) = lsqnonlin(f,30,-100,100,options);
            e(k)=eff(NTU(k),m_star_e_decoupled(k));
        end
        if is_dashed
            plot(NTU,e,'.r');
        else
            plot(NTU,e);
        end
        
        is_dashed = ~is_dashed;
    end  
end


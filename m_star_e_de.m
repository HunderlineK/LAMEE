function [ m_star_e ] = m_star_e_de( mse_in, mse_IWD, NTU, guess)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if mse_IWD==0
    m_star_e = mse_in;
else
    %m_star_e = mse_in*(1+1/2*mse_IWD*eff(NTU,guess)+1/6*(mse_IWD*eff(NTU,guess))^2);
    m_star_e = mse_in*(exp(mse_IWD*eff(NTU,guess))-1)/(mse_IWD*eff(NTU,guess));
end
end


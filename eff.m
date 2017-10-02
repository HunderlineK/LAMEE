function [ eff_ ] = eff( NTU, Cr )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if Cr == 1
        eff_ = NTU/(1+NTU);
    else
        eff_ = (1-exp(-NTU*(1-Cr)))/(1-Cr*exp(-NTU*(1-Cr)));
    end
    eff_ = max(0.001,eff_);
end


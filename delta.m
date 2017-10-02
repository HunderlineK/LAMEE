function [ d ] = delta(k, salt, design,Tair,Wair,Tsol,Xsol)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

 h_fg = 2.501 - 0.0024 * Tsol;
 cp_sol = HeatCapacity(salt,Tsol,Xsol);
 Wsol = W(VP(salt,Tsol,Xsol));
 
 dTair = design(1) * ( Tsol - Tair ) * 1/k;
 dWair = design(2) * ( Wsol - Wair ) * 1/k;
 dTsol = design(3) * ( dTair + h_fg * dWair );
 dXsol = -Xsol * design(3) * cp_sol * dWair * 1/1000;
 
 d(1) = dTair;
 d(2) = dWair;
 d(3) = dTsol;
 d(4) = dXsol;
end


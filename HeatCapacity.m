function HeatCapacity = HeatCapacity(desiccant, Z, C)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

Z = Z+273.15;

HCW = 88.7891 - 120.1958 * (Z / 228 - 1) ^ 0.02 - 16.9264 * (Z / 228 - 1) ^ 0.04 + 52.4654 * (Z / 228 - 1) ^ 0.06 + 0.10826 * (Z / 228 - 1) ^ 1.8 + 0.46988 * (Z / 228 - 1) ^ 8;
if strcmp(desiccant,'CaCl2')
    HeatCapacity = HCW * (1 - (1.63799 * C / 100 - 1.69002 * (C / 100) ^ 2 + 1.05124 * (C / 100) ^ 3) * ((58.225 * (Z / 228 - 1) ^ (0.02) - 105.6343 * (Z / 228 - 1) ^ (0.04) + 47.7948 * (Z / 228 - 1) ^ (0.06))));

elseif strcmp(desiccant,'LiCl')
    if C < 31 
        HeatCapacity = HCW * (1 - (1.4398 * (C / 100) - 1.24317 * (C / 100) ^ 2 - 0.1207 * (C / 100) ^ 3) * (58.5225 * (Z / 228 - 1) ^ 0.02 - 105.6343 * (Z / 228 - 1) ^ 0.04 + 47.7948 * (Z / 228 - 1) ^ 0.06));
    elseif C >= 31
        HeatCapacity = HCW * (1 - ((0.12825 + 0.62934 * (C / 100)) * (58.5225 * (Z / 228 - 1) ^ (0.02) - 105.6343 * (Z / 228 - 1) ^ (0.04) + 47.7948 * (Z / 228 - 1) ^ (0.06))));
    end
elseif strcmp(desiccant,'MgCl2')
    HeatCapacity = HCW + (-6304.3 + 3082.9 * C / 100 + 7.9 * (Z - 273.15) - 13.9 * 10 ^ (-3) * (Z - 273.15) ^ 2) * C / 100000;

elseif strcmp(desiccant,'LiBr')
    HeatCapacity = HCW + (-5277.4 + 2568.4 * C / 100 + 2.8 * (Z - 273.15) - 16.1 * 10 ^ (-3) * (Z - 273.15) ^ 2) * C / 100000;
else
    HeatCapacity = 0;
end

end


function W_ = W (VP)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global P_total

W_ = 0.62198 * (VP / (101325 / 1000- 0*VP)) * 1000;

end


function [ Twb ] = wetbulb( Tin, Win)

options = optimset('Display','off');
f = @(x) Tin - x - (W(VP('LiCl', x, 0.01))-Win)*2.4;
Twb = fsolve(f, 50,options);

end


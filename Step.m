function err = Step(k, salt, design, q, p)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

d = delta(k, salt, design, q(1),q(2),q(3),q(4));

err(1) = (q(1) - d(1) - p(1));
err(2) = (q(2) - d(2) - p(2));
err(3) = (q(3) - d(3) - p(3));
err(4) = (q(4) - d(4) - p(4));

end


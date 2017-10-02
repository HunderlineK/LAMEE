
function dy = dLAMEE(x,y)

global salt_ NTU_ NTUm_ Cr_

test = 'dLAMEE';

dy = zeros(4,1);
h_fg = min(2.501, abs(2.501 - 0.0024 * y(3)));
cp_sol = min(4.2, abs(HeatCapacity(salt_, y(3), y(4))));
Wsol = min(200, abs(W(VP(salt_, y(3), y(4)))));

dy(1) = -NTU_ * (y(3) - y(1));
dy(2) = -NTUm_ * (Wsol - y(2));
dy(3) = Cr_ * ( dy(1) + h_fg * dy(2));
dy(4) = -y(4) * Cr_ * cp_sol * dy(2) * 1/1000;

end


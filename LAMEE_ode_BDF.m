function out = LAMEE_ode (salt, NTU, NTUr, Cr, Tair, Wair, Tsol, Xsol)

global salt_ NTU_ NTUm_ Cr_ gTau gTair gWair gTsol gXsol

salt_ = salt;
NTU_ = NTU;
NTUm_ = NTU * NTUr;
Cr_ = Cr;

options = odeset('BDF','on');
[T,Y] = ode15s(@dLAMEE,[0,1],[Tair Wair Tsol Xsol],options);

subplot(1,4,1), plot(T,Y(:,1),'o')
subplot(1,4,2), plot(T,Y(:,2),'o')
subplot(1,4,3), plot(T,Y(:,3),'o')
subplot(1,4,4), plot(T,Y(:,4),'o')


out(1) = Y(length(Y),1);
out(2) = Y(length(Y),2);

gTau = T;
gTair = Y(:,1);
gWair = Y(:,2);
gTsol = Y(:,3);
gXsol = Y(:,4);
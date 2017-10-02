function out = outlet (salt, NTU, NTUr, Cr, Tair_in, Wair_in, Tsol_in, Xsol_in)

NTUm = NTU * NTUr;

design(1) = NTU;
design(2) = NTUm;
design(3) = Cr;

k = 100;

Tair = 1:k;
Wair = 1:k;
Tsol = 1:k;
Xsol = 1:k;

Tair(1) = Tair_in;
Wair(1) = Wair_in;
Tsol(1) = Tsol_in;
Xsol(1) = Xsol_in;

for i=1:k
 options = optimset('Display','off');
 s = @(x)Step(k, salt, design, x,[Tair(i),Wair(i),Tsol(i),Xsol(i)]);
 x = fsolve(@(x) s(x),[Tair(i),Wair(i),Tsol(i),Xsol(i)],options);
 
 Tair(i+1) = x(1);
 Wair(i+1) = x(2);
 Tsol(i+1) = x(3);
 Xsol(i+1) = x(4);

 
end

subplot(1,4,1), plot(0:k,Tair,'.')
subplot(1,4,2), plot(0:k,Wair,'.')
subplot(1,4,3), plot(0:k,Tsol,'.')
subplot(1,4,4), plot(0:k,Xsol,'.')

Tair = Tair(k+1);
Wair = Wair(k+1);
Tsol = Tsol(k+1);
Xsol = Xsol(k+1);

out(1) = Tsol;
out(2) = Xsol;
%out(3) = Tair;
%out(4) = Wair;


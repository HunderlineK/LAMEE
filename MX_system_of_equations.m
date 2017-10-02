% velocity profile v(z'=z/z0)= a z'^2 + b z' + c
% v(z'=0=0/z0) = 0 -> c = 0
% v(z'=2=2*z0/z0) = 0 -> 4 a + 2 b = 0 -> 2 a + b = 0
% int(v(z'=z/z0),from 0 to z0) = v_ave -> 1/3*a+b/2 = v_ave 
%-> 2 a + 3 b =6 v_ave -> 2 b = 6 v_ave -> b = 3 v_ave -> a = -3/2 v_ave
%-> v(z')/v_ave = -3/2 z'^2 + 3 z'

clc;
clear all

n = 500;
m = 500;

v = @(z) -3/2.*(z./n-0.5/n).^2 + 3.*(z./n-0.5/n); %velocity profile

z0 = 0.0025;
x0 = 0.5;
SA = 0.1; %surface area
dA = SA/m;

mair = 0.00148;
Cr = 0.5;

Tsol_in = 24;
Wair_in = 17.5;
Xsol_in = 32;

Dair = 2.52*10^-5;
rho = 1.15;
cp_air = 1006;
h_fg = 2501-2.4*Tsol_in;

Cr_m = Cr*h_fg/cp_air; %Cr for moisture transfer

dWT = (W(VP('LiCl', Tsol_in, Xsol_in))-W(VP('LiCl', Tsol_in-1, Xsol_in)));
Wsol_in = W(VP('LiCl', Tsol_in, Xsol_in));

Um_sa = rho * Dair/((z0/n)/2); %kair/(z0/(n/2)); %U between solution and adjacent air
Um_aa_z =  rho * Dair/(z0/n); %U between adjacent air strips oriented in the x direction

A = spalloc( (1+n) * m, (1+n) * m, 4 * (1+n) * m);
%A = zeros( m*(n+1), m*(n+1));
b = zeros( (1+n) * m, 1);

%solution inlet condition
A( 1, 1) = dWT*Um_sa*dA/mair + 1/Cr_m;
A( 1, 2*m) = -Um_sa*dA/mair;
b(1) = Tsol_in * (1/Cr_m + dWT*Um_sa*dA/mair) - Wsol_in*Um_sa*dA/mair;

%/air inlet condition: start
%/nodes adjacent to solution
A( m+1, m) = dWT*Um_sa*dA/mair;
if n>=2
    A( m+1, 2*m+1) = Um_aa_z*dA/mair;
    A( m+1, m+1) = -Um_sa*dA/mair -Um_aa_z*dA/mair - 1/n*v(1);
    b( m+1) = Tsol_in*dWT*Um_sa*dA/mair -Um_sa*dA/mair*Wsol_in -Wair_in*v(1)/n;
else
    A( m+1, m+1) = -Um_sa*dA/mair - 1/n;
    b( m+1) = Tsol_in*dWT*Um_sa*dA/mair -Um_sa*dA/mair*Wsol_in -Wair_in;
end
for i = 2:m
    A( m+i, m - i + 1) = dWT*Um_sa*dA/mair;
    if n>=2
        A( m+i, 2*m+i) = Um_aa_z*dA/mair;
        A( m+i, m+i) = -Um_sa*dA/mair - Um_aa_z*dA/mair - 1/n*v(1);
        A( m+i, m+i-1) = 1/n*v(1);
    else
        A( m+i, m+i) = -Um_sa*dA/mair - 1/n;
        A( m+i, m+i-1) = 1/n;
    end
    b( m+i) = Tsol_in*dWT*Um_sa*dA/mair -Wsol_in*Um_sa*dA/mair;
end
%nodes adjacent to solution/
%/nodes adjacent to adiabatic wall
if n>=2
    A( n*m+1, (n-1)*m+1) = Um_aa_z*dA/mair;
    A( n*m+1, n*m+1) = - Um_aa_z*dA/mair - 1/n*v(n);
    b( n*m+1) = -Wair_in/n;
    for i = 2:m
        A( n*m+i, (n-1)*m + i) = Um_aa_z*dA/mair;
        A( n*m+i, n*m+i) = - Um_aa_z*dA/mair - 1/n*v(n);
        A( n*m+i, n*m+i-1) = 1/n*v(n);
    end   
end
%nodes adjacent to adiabatic wall/

%/nodes between
for j = 2:n-1   
    A( j*m+1, (j-1)*m+1) = Um_aa_z*dA/mair;
    A( j*m+1, j*m+1) = -Um_aa_z*dA/mair -Um_aa_z*dA/mair - 1/n*v(j);
    A( j*m+1, (j+1)*m+1) = Um_aa_z*dA/mair;
    b( j*m+1) = -Wair_in/n*v(j);
    for i = 2:m
        A( j*m+i, (j-1)*m+i) = Um_aa_z*dA/mair;
        A( j*m+i, j*m+i) = -Um_aa_z*dA/mair -Um_aa_z*dA/mair - 1/n*v(j);
        A( j*m+i, (j+1)*m+i) = Um_aa_z*dA/mair;
        A( j*m+i, j*m+i-1) = 1/n*v(j);    
    end
end
%nodes between/
%air inlet condition: finish/
    
%solution channel nodes
for i = 2:m
    A( i, i) = dWT*Um_sa*dA/mair + 1/Cr_m;
    A( i, 2*m+1-i) = -Um_sa*dA/mair;
    A( i, i-1) = -1/Cr_m;
    b(i) = Tsol_in * dWT*Um_sa*dA/mair - Wsol_in*Um_sa*dA/mair;
end

A=sparse(A);
b=sparse(b);

out= A\b;

x = (1:m)/m;

Tsol = out(1:m);
for i=1:m
    Wsol(i) = Wsol_in + dWT*(Tsol(i)-Tsol_in);
end

for j=1:n
    Wair(j,:) = out(j*m+1:(j+1)*m);
end

for i=1:m
    Wair_average(i) = sum(Wair(:,i).*v(1:n)')/n;
end


hold off;
plot((z0/n)*(1:n), Wair_average(:));
hold on;
plot((z0/n)*(n:-1:1), Wsol(:));

figure
surf(Wair,'EdgeColor','none');
view(2);
colorbar






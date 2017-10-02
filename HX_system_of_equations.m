% velocity profile v(z'=z/z0)= a z'^2 + b z' + c
% v(z'=0=0/z0) = 0 -> c = 0
% v(z'=2=2*z0/z0) = 0 -> 4 a + 2 b = 0 -> 2 a + b = 0
% int(v(z'=z/z0),from 0 to z0) = v_ave -> 1/3*a+b/2 = v_ave 
%-> 2 a + 3 b =6 v_ave -> 2 b = 6 v_ave -> b = 3 v_ave -> a = -3/2 v_ave
%-> v(z')/v_ave = -3/2 z'^2 + 3 z'


clc;
clear all

n = 1000;
m = 1000;

kair = 0.0257;
z0 = 0.0025;
x0 = 0.5;
SA = 0.1; %surface area
dA = SA/m;

Cair = 1.5;
Cr = 1;
Tsol_in = 24;
Tair_in = 34;


U_sa = kair/((z0/n)); %kair/(z0/(n/2)); %U between solution and adjacent air
U_aa_z = kair/(z0/n); %U between adjacent air nodes

v = @(z) -3/2.*(z./n-0.5/n).^2 + 3.*(z./n-0.5/n); %velocity profile

A = spalloc( (1+n) * m, (1+n) * m, 4 * (1+n) * m);
b = zeros( (1+n) * m, 1);

%solution inlet condition
A( 1, 1) = U_sa*dA/Cair + 1/Cr;
A( 1, 2*m) = -U_sa*dA/Cair;
b( 1) = Tsol_in * 1/Cr;

%/air inlet condition: start
%/nodes adjacent to solution
A( m+1, m) = U_sa*dA/Cair;
if n>=2
    A( m+1, 2*m+1) = U_aa_z*dA/Cair;
    A( m+1, m+1) = -U_sa*dA/Cair -U_aa_z*dA/Cair - 1/n*v(1);
else
    A( m+1, m+1) = -U_sa*dA/Cair - 1/n*v(1);
end
b( m+1) = -Tair_in/n*v(1);

for i = 2:m
    A( m+i, m - i + 1) = U_sa*dA/Cair;
    if n>=2
        A( m+i, 2*m+i) = U_aa_z*dA/Cair;
        A( m+i, m+i) = -U_sa*dA/Cair - U_aa_z*dA/Cair - 1/n*v(1);
    else
        A( m+i, m+i) = -U_sa*dA/Cair - 1/n*v(1);
    end
    A( m+i, m+i-1) = 1/n*v(1);
end
%nodes adjacent to solution/
%/nodes adjacent to adiabatic wall
if n>=2
    A( n*m+1, (n-1)*m+1) = U_aa_z*dA/Cair;
    A( n*m+1, n*m+1) = - U_aa_z*dA/Cair - 1/n*v(n);
    b( n*m+1) = -Tair_in/n;
    for i = 2:m
        A( n*m+i, (n-1)*m + i) = U_aa_z*dA/Cair;
        A( n*m+i, n*m+i) = - U_aa_z*dA/Cair - 1/n*v(n);
        A( n*m+i, n*m+i-1) = 1/n*v(n);
    end   
end
%nodes adjacent to adiabatic wall/

%/nodes between
for j = 2:n-1   
    A( j*m+1, (j-1)*m+1) = U_aa_z*dA/Cair;
    A( j*m+1, j*m+1) = -U_aa_z*dA/Cair -U_aa_z*dA/Cair - 1/n*v(j);
    A( j*m+1, (j+1)*m+1) = U_aa_z*dA/Cair;
    b( j*m+1) = -Tair_in/n*v(j);
    for i = 2:m
        A( j*m+i, (j-1)*m+i) = U_aa_z*dA/Cair;
        A( j*m+i, j*m+i) = -U_aa_z*dA/Cair -U_aa_z*dA/Cair - 1/n*v(j);
        A( j*m+i, (j+1)*m+i) = U_aa_z*dA/Cair;
        A( j*m+i, j*m+i-1) = 1/n*v(j);    
    end
end
%nodes between/
%air inlet condition: finish/
    
%solution channel nodes
for i = 2:m
    A( i, i) = U_sa*dA/Cair + 1/Cr;
    A( i, 2*m+1-i) = -U_sa*dA/Cair;
    A( i, i-1) = -1/Cr;
end

A=sparse(A);
b=sparse(b);

T = A\b;

x = (1:m)/m;

Tsol = T(1:m);

for j=1:n
    Tair(j,:) = T(j*m+1:(j+1)*m);
end

for i=1:m
    Tair_average(i) = sum(Tair(:,i).*v(1:n)')/n;
end;

Tair_average(i)

hold off;
plot((z0/n)*(1:n), Tair_average(:));
hold on;
plot((z0/n)*(n:-1:1), Tsol(:));

figure
surf(Tair,'EdgeColor','none');
view(2);
colorbar

 
        
%psychometric
m_ = 0;
n_ = 60;
Wsat_T = ones(n_-m_+1,1);
T = ones(n_-m_+1,1);
for k = 0:60
    T(k+1) = k;
    Wsat_T(k+1) = W(VP('LiCl', T(k+1), 0.01));
end
hold on
plot(T,Wsat_T,'r');
hold off
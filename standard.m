%STD results
for k=m:n
    es_std(k) = eff(NTU_(k),Cr_(k));
    dTair_std(k) = es_std(k) * (Tsol_in_(k)-Tair_in_(k));
    
    m_star(k) = cp_sol(k)* Cr_(k) / 2.4;
    em_std(k) = eff(NTU_m(k),m_star(k));
    dWair_std(k) = em_std(k) * (Wsol_in_(k)-Wair_in_(k));
end
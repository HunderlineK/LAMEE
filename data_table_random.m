%Creates random data points and invokes the numerical solution to calculate the numerical outlet data and the simplified-extended effectiveness models to calculate the estimated outlet data, for each data point for the stead state operating conditions

clc
% clears the MATLAB command screen
clear
% clears all the variables

global gTau gTair gWair gTsol gXsol;
% defines the global variables that will be accessible to all the functions
% gTau is the variable that holds the dimensionless points along the exchanger at which discretization is applied
% gTair, gWair, gTsol and gXsol are, respectively, the temperature and humidity ratio profiles of the air stream and temeprature and concetration profiles of the solution stream, along the exchange

m = 1
n = 5000
% total number of data points equals to n-m+1

for p = m:n
% the loop to generate the random inlet operating conditions and design paramerters for the data points
   Cr_(p)= random('unif',0.1,0.5);
   % assigns the heat capacity rate ratio
   NTU_(p) = random('unif',1,6);
   % assigns the NTU
   NTU_ratio_(p) = random('unif',0.2,0.8);
   % assigns the ratio NTU_m/NTU
   NTU_m(p) = NTU_(p) * NTU_ratio_(p);
   % calculates and assigns NTU_m
   
   if random('unif',0,2) >= 1
       salt_{p} = 'LiCl';
   else
       salt_{p} = 'MgCl2';
   end
   % assigns the salt type
   
   Tsol_in_(p) = random('unif',15,25);
   % assigns the inlet solution temperature
   Xsol_in_(p) = random('unif',25,35);
   % assigns the inlet solution concentration
   Wsol_in_(p) = W(VP(salt_{p}, Tsol_in_(p), Xsol_in_(p)));
   % calculates and assigns the inlet solution equilibirum humidity ratio
   
   Tair_in_(p) = random('unif', 30,40);
   % assigns the inlet air temperature
   Wair_in_(p) = random('unif',12,24);
   % assigns the inlet air humidity ratio
   
   h_fg(p) = 2.501 - 0.0024 * Tsol_in_(p);
   % calculates and assigns the specific heat of evaporation
   cp_sol(p) = HeatCapacity(salt_{p}, (Tsol_in_(p)), Xsol_in_(p));
   % calculates and assigns the specific heat of solution
   h_star(p) = h_fg(p) * ( Wsol_in_(p) - Wair_in_(p))/(Tsol_in_(p) - Tair_in_(p));
   % calculates and assigns h_star

end

for p=m:n
% the loop for numerical calculations and also estimations of the simplified-extended
% method for outlet air operating conditions

   % start of the calculations for the estimations provided by the simplified-extended effectiveness models
   mse_in_MX(p) = 0.058 * Wsol_in_(p) * Cr_(p) * h_fg(p);
   % calculates and assigns mse_in_MX
   mse_IWD_MX(p) = 0.058 * (Wair_in_(p)-Wsol_in_(p)) * Cr_(p) * h_fg(p);
   % calculates and assigns mse_IWD_MX
   f = @(x) (x - m_star_e_de(mse_in_MX(p), mse_IWD_MX(p), NTU_m(p),x));
   % assigns the function f which will be solved to find x, where x is the m_star_e_MX
   fval = 0.01;
   % sets the initial residue of the solver
   trial = 1;
   % initializes the number of trials to 1
   options = optimset('Display','off');
   % sets the display option to off
   while fval>=0.01 && trial <= 150
   % the loop that will continue for at most 150 trials until fval is less than 0.01
       if trial == 1
       % for the first trial, the starting point of the solver is set to x=20
       % m_star_e_decoupled is the variable that represents the solution for x, which in turn is m_star_e_MX
           [m_star_e_decoupled(p),fval] = lsqnonlin(f, 20,-100,100,options);
       else
       % for the other trials, the starting point of the solver is randomly selected between -100 and 100
           [m_star_e_decoupled(p),fval] = lsqnonlin(f, random('unif',-100,100),-1e3,1e3,options);
       end
       trial = trial+1;
       % increases the trial number after each attempt to solve for x
   end
   
   eff_r(p) = eff(NTU_m(p),0)/eff(NTU_(p),0);
   % calculates and assigns the approximation for the ratio e_m/e_s

   
   Cr_e_emp(p) = Cr_(p) * (1 + h_star(p) * eff_r(p));
   % calculates and assigns the estimated Cr_e
   eff_s_emp(p) = eff(NTU_(p), Cr_e_emp(p));
   % calculates and assigns eff_s estimated based on the estimated Cr_e
   dTair_emp(p) = eff_s_emp(p)*(Tsol_in_(p) - Tair_in_(p));
   % calculates and assigns dTair (change in the temperature of the air stream) estimated based on the estimated eff_s
   
   
   mse_in(p) = 0.058 * Wsol_in_(p) * Cr_(p) * h_fg(p) * (1+1/(h_star(p) * eff_r(p)));
   % calculates and assigns mse_in
   mse_IWD(p) = 0.058 * (Wair_in_(p)-Wsol_in_(p)) * Cr_(p) * h_fg(p) * (1 + 1 /(h_star(p) * eff_r(p)));
   % calculate and assigns mse_IWD   
   f = @(x) (x - m_star_e_de(mse_in(p), mse_IWD(p), NTU_m(p),x));
   % assigns the function f which will be solved to find x, where x is the m_star_e
   
   fval = 0.1;
   trial = 1;
   options = optimset('Display','off');
   while fval>=0.01 && trial <= 150
       if trial == 1
       % m_star_e_coupled is the variable that represents the solution for x, which in turn is m_star_e
           [m_star_e_coupled(p),fval] = lsqnonlin(f, 20,-100,100,options);
       else
           [m_star_e_coupled(p),fval] = lsqnonlin(f, random('unif',-100,100),-1e3,1e3,options);
       end
       trial = trial+1;
   end
   if fval>0.1
   % if the residua of the solver is greater than 0.1, sets fail flag to 1
       failed_emp(p)=1;
   else
       failed_emp(p)=0;
   end

   eff_m_emp(p) = eff(NTU_m(p), m_star_e_coupled(p));
   % calculates and assigns eff_m based on the estimated values for m_star_e
   dWair_emp(p) = eff_m_emp(p)*(Wsol_in_(p) - Wair_in_(p));
   % calculates and assigns dWair (change in the humidity ratio of the air stream) based on the estimated eff_m
   dHair_emp(p) = dTair_emp(p) + h_fg(p) * dWair_emp(p);
   % calculates and assigns dHair based on the estimated  dTair and dWair

   Tair_out_emp (p) = Tair_in_(p) + dTair_emp(p);
   % calcualtes and assigns Tair_out based on the estimated dTair
   Wair_out_emp (p) = Wair_in_(p) + dWair_emp(p);
   % calculates and assigns Wair_out based on the estimated dWair
   % the calculations for the estimations provided by the simplified-extended effectiveness models end here
   
   
   % Calculations for the numerical model start here
   fval = 0.1;
   trial = 1;
   while fval>=0.1 && trial <= 50
       if trial == 1
           % 1st guess solution = empirical output
           x01 = Tair_out_emp (p);
           x02 = Wair_out_emp (p);
           trial_emp = [p, x01, x02]
       else
           % other guess solutions = randomly selected
           x01 = max(0.5,random('norm',Tair_out_emp (p),15));
           Wsat_x01 = W(VP(salt_{p}, x01, 0.1));
           x02 = max(0.5,random('norm',Wair_out_emp (p),15));
           trial_ran = [p,trial, x01, x02]
       end
       [x,fval] = lsqnonlin(@(x) LAMEE_ode(salt_{p}, NTU_(p) , NTU_ratio_(p), ...
           Cr_(p),...
           x(1),x(2), Tsol_in_(p), Xsol_in_(p))... 
           -[ Tair_in_(p), Wair_in_(p)],[x01 ,  x02],[0.5 0.5], [80 W(VP(salt_{p}, 80, 0.1))],options);
       % solves for x(1) and x(2), which are Tair_out_(p) and Wair_out_(p). The solution for  x(1) and x(2)
       % results in Tair_in and Wair_in that are the same as (or closest to) Tair_in_(p) and Wair_in_(p)
       trial = trial+1;
   end
   if fval>1
       failed_num(p)=1
   else
       failed_num(p)=0;
   end
   
   Tair_out(p) = real(gTair(1));
   % the variable gTair holds the temperature profile of the air along the exchanger
   Wair_out(p) = real(gWair(1));
   % the variable gWair holds the humidity ratio profile of the air along the exchanger
   
   dTair_num(p) = Tair_out(p) - Tair_in_(p);
   dWair_num(p) = Wair_out(p) - Wair_in_(p);
   
   Table{p,1} = gTau;
   Table{p,2} = gTair;
   Table{p,3} = gWair;
   Table{p,4} = gTsol;
   Table{p,5} = gXsol;
   % store profiles of Tair, Wair, Tsol and Xsol, and gTau (the points along the exchangers which data is available)
  
end
   

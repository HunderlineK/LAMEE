%% PIREG PI-regulator. Ger regulatorv�rde.
%%
%% ut =  pireg(ipart,meas,set,k,itime)

function ut =  pireg(ipart,meas,set,k,itime)

e = set - meas;

ut = k*(e + 1/itime*ipart);


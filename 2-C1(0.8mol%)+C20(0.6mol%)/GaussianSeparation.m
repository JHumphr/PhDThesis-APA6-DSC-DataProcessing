function k = GaussianSeparation(x,y,Hc)
%%%This function takes in x, y of time and baseline heatflow and returns
fun = @(k)SAGauss(x,k);

 k0 = [0.23 11 1 1 0.3 12 1 1];  
% 
 Kmin= [0 10.5 0.85 0.75 0 14 0.5 0.5];
 Kmax= [0.6 12 inf inf inf 18 inf inf];

k = fmincon(@(k)y-fun,k0,[],[],[],[],Kmin,Kmax,@CONSTRAIN); 

end


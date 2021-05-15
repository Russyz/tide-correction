function ng = Atm_refractive_index(wavelength,baroPressure,temperature,RH)
wavelength = 1064;
T=temperature+273.15;
Ps=6.11*10^(((7.5*(T-273.15)))/((237.3+(T-273.15))));
e=Ps*RH/100.0;
ng = 1+1e-6*(273.15*baroPressure/(1013.25*T)*(287.6155+4.88660/wavelength^2+0.06800/wavelength^4)-11.27*e/T);
end

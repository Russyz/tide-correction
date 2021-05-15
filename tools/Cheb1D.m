%--------------------------------------------------------------------------
%
% Chebyshev approximation of 1-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     C      coeff
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function ChebApp = Cheb1D(t, N, Ta, Tb, C)
% Check validity
if ( (t<Ta) || (Tb<t) )
    error('ERROR: Time out of range in Cheb3D::Value\n');
end
% Clenshaw algorithm
tau = (2*t-Ta-Tb)/(Tb-Ta);  
f1 = zeros(1,1);
f2 = zeros(1,1);
for i=N:-1:2
    old_f1 = f1;
    f1 = 2*tau*f1-f2+C(i);
    f2 = old_f1;
end
ChebApp = tau*f1-f2+C(1);
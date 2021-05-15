function [yi, ypi, yppi] = lagint(x, y, xi, c)
% LAGINT 1-D piecewise lagrange interpolation
%    LAGINT(X,Y,XI,C) interpolates to find YI, the values of the underlying
%    function Y at the points in the array XI, using piecewise lagrange
%    interpolation.  X and Y must be vectors of length N.
%
%    C specifies the number of data points to use in the interpolation. 
%    The default is to use all points.
%
%    [YI,YPI,YPPI] = LAGINT() also returns the interpolated first and
%    second derivatives of the underlying function Y at points XI.
% Joe Henning - Fall 2011
if (nargin < 4)
   c = 0;
end
n = length(x);
if n ~= length(y)
   fprintf('??? Bad input to lagint ==> X and Y must be of the same length\n');
   yi = [];
   ypi = [];
   yppi = [];
   return;
end
if (n < c)
   fprintf('??? Bad input to lagint ==> C <= N\n');
   yi = [];
   ypi = [];
   yppi = [];
   return
end
for i = 1:length(xi)
   % Find the right place in the table by means of a bisection.
   klo = 1;
   khi = n;
   while (khi-klo > 1)
      k = fix((khi+klo)/2.0);
      if (x(k) > xi(i))
         khi = k;
      else
         klo = k;
      end
   end
   
   h = x(khi) - x(klo);
   if (h == 0.0)
      fprintf('??? Bad input to lagint ==> X values must be distinct\n');
      yi(i) = NaN;
      ypi(i) = NaN;
      yppi(i) = NaN;
      continue;
   end
   % Evaluate lagrange polynomial
   yi(i) = 0;
   ypi(i) = 0;
   yppi(i) = 0;
   if (c == 0)
      for k = 1:n
         term = y(k);
         termp = 0;
         termp2 = 0;
         for m = 1:n
            if (k ~= m)
               term = term*(xi(i)-x(m))/(x(k)-x(m));
               prod = 1;
               for j = 1:n
                  if ((k ~= j) && (m ~= j))
                     prod = prod*(xi(i)-x(j))/(x(k)-x(j));
                  end
               end
               termp = termp + y(k)*prod/(x(k)-x(m));
               termp3 = 0;
               for j = 1:n
                  if ((k ~= j) && (m ~= j))
                     prod = 1;
                     for l = 1:n
                        if ((k ~= l) && (m ~= l) && (j ~= l))
                           prod = prod*(xi(i)-x(l))/(x(k)-x(l));
                        end
                     end
                     termp3 = termp3 + prod/(x(k)-x(j));
                  end
               end
               termp2 = termp2 + y(k)*termp3/(x(k)-x(m));
            end
         end
         yi(i) = yi(i) + term;
         ypi(i) = ypi(i) + termp;
         yppi(i) = yppi(i) + termp2;
      end
   else
      if (mod(c,2) == 0)   % even
         c2 = c/2;
         if (klo < c2)
            klo = c2;
         end
         if (klo > n-c2)
            klo = n-c2;
         end
         khi = klo + 1;
         for k = klo-(c2-1):klo+c2
            term = y(k);
            termp = 0;
            termp2 = 0;
            for m = klo-(c2-1):klo+c2
               if (k ~= m)
                  term = term*(xi(i)-x(m))/(x(k)-x(m));
                  prod = 1;
                  for j = klo-(c2-1):klo+c2
                     if ((k ~= j) && (m ~= j))
                        prod = prod*(xi(i)-x(j))/(x(k)-x(j));
                     end
                  end
                  termp = termp + y(k)*prod/(x(k)-x(m));
                  termp3 = 0;
                  for j = klo-(c2-1):klo+c2
                     if ((k ~= j) && (m ~= j))
                        prod = 1;
                        for l = klo-(c2-1):klo+c2
                           if ((k ~= l) && (m ~= l) && (j ~= l))
                              prod = prod*(xi(i)-x(l))/(x(k)-x(l));
                           end
                        end
                        termp3 = termp3 + prod/(x(k)-x(j));
                     end
                  end
                  termp2 = termp2 + y(k)*termp3/(x(k)-x(m));
               end
            end
            yi(i) = yi(i) + term;
            ypi(i) = ypi(i) + termp;
            yppi(i) = yppi(i) + termp2;
         end
      else   % odd
         c2 = floor(c/2);
         if (klo < c2+1)
            klo = c2+1;
         end
         if (klo > n-c2)
            klo = n-c2;
         end
         khi = klo + 1;
         for k = klo-c2:klo+c2
            term = y(k);
            termp = 0;
            termp2 = 0;
            for m = klo-c2:klo+c2
               if (k ~= m)
                  term = term*(xi(i)-x(m))/(x(k)-x(m));
                  prod = 1;
                  for j = klo-c2:klo+c2
                     if ((k ~= j) && (m ~= j))
                        prod = prod*(xi(i)-x(j))/(x(k)-x(j));
                     end
                  end
                  termp = termp + y(k)*prod/(x(k)-x(m));
                  termp3 = 0;
                  for j = klo-c2:klo+c2
                     if ((k ~= j) && (m ~= j))
                        prod = 1;
                        for l = klo-c2:klo+c2
                           if ((k ~= l) && (m ~= l) && (j ~= l))
                              prod = prod*(xi(i)-x(l))/(x(k)-x(l));
                           end
                        end
                        termp3 = termp3 + prod/(x(k)-x(j));
                     end
                  end
                  termp2 = termp2 + y(k)*termp3/(x(k)-x(m));
               end
            end
            yi(i) = yi(i) + term;
            ypi(i) = ypi(i) + termp;
            yppi(i) = yppi(i) + termp2;
         end
      end
   end
end
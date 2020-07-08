function Sy = makeSimpVec(dy,Ny)
% makeSimpVec makes a vector Sy that can be used for integration by
% Simpson's rule for a function with step dy and Ny entries.
% The error of the integration is on the order dy^5*f''''(y).
% Usage: Sy = makeSimpVec(dy,Ny);
%        integral = sum(f.*Sy);
if(mod(Ny,2)==1 && Ny>=3)
    Sy = (dy/3)*[1 repmat([4 2], [1 (Ny-3)/2]) 4 1];
elseif(mod(Ny,3)==1 && Ny>=4)
    Sy = (dy/8)*[3 repmat([9 9 6], [1 (Ny-4)/3]) 9 9 3];
elseif(Ny>=4)
    Sy = dy*[(1/3)*[1 repmat([4 2], [1 (Ny-4)/2]) 4] 1/3+1/2 1/2];
else
    Sy = (dy/2)*ones(1,Ny);
end
% All were tested and give accurate results!

%error = (1/90)h^5*f(4)
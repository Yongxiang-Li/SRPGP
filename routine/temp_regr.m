function [ H ] = temp_regr( temp, func, x )
% new basis: Bspline + temperature
    
    B = func(x);
    T = [temp temp.^2];
    H = [T B'];

end


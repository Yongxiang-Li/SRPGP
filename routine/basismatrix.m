function B = basismatrix(d, p, knots, U) 
%  
% Function Name: 
%  
%   bspeval - Evaluate a univariate B-Spline. 
%  
% Calling Sequence: 
%  
%   p = bspeval(d,c,k,u) 
%  
% Parameters: 
%  
%   d	: Degree of the B-Spline. 
%  
%   knots	: Knot sequence, row vector of size nk. 
%  
%   U	: Parametric evaluation points, row vector of size nu. 
%  
%   B	: Evaluated points, matrix of size (dim,nu) 
%  
% Description: 
%  
%   Evaluate a univariate B-Spline.
nu = numel(U); 
B = sparse(p, nu); 
for col=1:nu 
    s = findspan(p-1, d, U(col), knots);
    N = basisfun(s,U(col),d,knots);
    B(s-d+1:s+1, col) = N;
end


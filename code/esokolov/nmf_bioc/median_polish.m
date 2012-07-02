function [x, grand_effect, row_effect, column_effect, resid] = median_polish(x, tol, maxiter)
% MEDIAN_POLISH Fits an additive model using Tukey's median polish procedure. 
%   MEDIAN_POLISH(X) fits an additive model to a 2D input matrix, X, assuming a 
%   decomposition of (constant + rows + columns).
%
%   Median polishing is useful for removing spatial trends in the data by 
%   alternately removing medians from the rows and columns of the data. The 
%   algorithm proceeds as follows. Considering the rows first, for each row the 
%   row median is subtracted from every element in that row. Then, for each 
%   column, the median of the revised numbers is then subtracted from every 
%   element in that column. This process is repeated until all convergence, or 
%   a maximum number of iterations.
%
%   In this implementation, NaNs are assumed to be missing data.
%
%   Syntax:
%
%   p = median_polish(X);
%   p = median_polish(X, tol);
%   p = median_polish(X, tol, maxiter);
%   [p, ge] = median_polish(...);
%   [p, ge, re] = median_polish(...);
%   [p, ge, re, ce] = median_polish(...);
%   [p, ge, re, ce, resid] = median_polish(...);
%   
%   X       : A 2D input matrix
%   tol     : Tolerance threshold for stopping iterations
%   maxiter : The maximum number of iterations
%
%   p       : The median polished matrix.
%   ge      : The fitted constant term (also known as the Grand Effect).
%   re      : The fitted row effect.
%   ce      : The fitted column effect.
%   resid   : Residuals
%
%   Example: If X = [13 17 26 18 29;
%                    42 48 57 41 59;
%                    34 31 36 22 41];
%
%   then median_polish(X) returns [ 0 -1  0  7  0;
%                                  -1  0  1  0  0;
%                                   9  1 -2 -1  0]
%
%   Class support for input X:
%      float: double, single
%
%   Reference: Tukey, John W. Exploratory Data Analysis, Addison-Wesley 1977 
%              ISBN 0-201-07616-0
%
%   See also MEDIAN, NANMEDIAN.

%   Adam Auton 2009
%   $Revision: 1.0 $  $Date: 2009/03/28 14:00:00 $
%   Tested on MATLAB Release R2008a

Narg = nargin;
if (Narg < 2)
    tol = 0.01;
end
if (Narg < 3)
    maxiter = 100;
end

S = size(x);
N_dim = length(S);

% Parameter Check
if (N_dim ~= 2)
    error('Require a 2 dimensional matrix as input');
end
if (tol <= 0)
    error('Convergence Tolerance (tol) must be positive');
end
if (maxiter < 1)
    error('Maximum number of iterations must greater or equal to 1');
end

xorig = x;
grand_effect = NaN;
row_effect = zeros(1, S(2));
column_effect = zeros(S(1), 1);

for it=1:maxiter    
    
    % Row Sweep
    rm = nanmedian(x, 1);
    x = bsxfun(@minus, x, rm);

    row_effect = (row_effect + rm);
    rmm = nanmedian(row_effect);

    
    % Column Sweep
    cm = nanmedian(x, 2);
    x = bsxfun(@minus, x, cm);
    
    column_effect = column_effect + cm;
    cmm = nanmedian(column_effect);

    % Grand Effect
    ge_prime = grand_effect;
    grand_effect = rmm + cmm;
    
    d = abs(ge_prime - grand_effect);
    if (d < tol)
        break;
    end
end

resid = xorig - x;  % Calculate Residuals

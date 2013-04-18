% routine to find initial distribution of points with latin hypercube sampling
% ninit is number of points to include
% amin is vector of lower bounds for params
% amax is vector of upper bounds for params
% N is number of params

function init_pts = initial_pts_Nd(ninit, amin, amax, N)

Ssc = lhsamp(ninit, N);
% scale points to bounds for parameters
for i = 1:N
    S(:,i) = Ssc(:,i).*(amax(i) - amin(i)) + amin(i);
end

init_pts = S;
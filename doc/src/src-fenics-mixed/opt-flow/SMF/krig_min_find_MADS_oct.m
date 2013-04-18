% Alison Marsden
% function krig_min_find_MADS
% uses DACE package to find minimum of kriging fit for N dimensions
% returns minimum surrogate points (next_pts) and maximum mse point
% (max_mse_pt)
% use MADS polling to find additional search point
% last update June 28, 2006

function [next_pts, min_est, max_mse_pt] = krig_min_find_MADS_oct(S,Y,curr_bestx,theta,upb,lob,N,amin,amax, spc, delta)

% call dacefit to do kriging fit and return model
[dmodel, perf] = dacefit(S,Y, @regpoly0, @corrgauss, theta, lob, upb);
disp('dacefit called')
clear perf theta lob upb Y S

% find search point based on polling surrogate around current best point,
% taking lowest value
poll_pts = MADS_poll_ptsNd3_oct(curr_bestx,N,delta,spc,amin,amax)
%poll_pts = find_poll_ptsNd(curr_bestx,N,spc,amin,amax);
if isempty(poll_pts) == 0
    for i = 1:size(poll_pts,1)
        [Jtest(i)] = predictor(poll_pts(i,:), dmodel)
    end
    [low_Jtest,i] = min(Jtest)
    poll_xpt = poll_pts
end
% find  surrogate minimum points
% several surrogate minimum points are found, and sorted in order of lowest to highest cost
% function value 
% number of surrogate minimum points = nsm
nsm = 1;
strfun = 'predictor';
xlow = [amin]';
xup = [amax]';
for n = 1:nsm
    options = es_options_new('TolX',1e-4,'SigmaFacStart',1,'LBound',amin','UBound',amax');
    es = es_new(N, xlow, xup,options);     % generate an ES struct
    es = es_run_mod(es, strfun, inf, dmodel);       % minimize strfun with ES
    result = es_get(es, 'result'); % get result
    xevotmp(:,n) = result.x;               % return current best point
    [fevotmp(n)] = predictor(xevotmp(:,n), dmodel);
end
[fevos,I] = sort(fevotmp);
next_pts = zeros(1,N);
for j = 1:nsm
    fevo = fevos(j);
    nn = I(j);
    xevo = xevotmp(:,nn);
    for i=1:N
        interv(i) = (amax(i)-amin(i))/(spc(i)-1);
    end
    % find nearest mesh point
    for n=1:N
        ind1 = 1+floor((xevo(n)-amin(n))/interv(n));
        ind2 = 1+ceil((xevo(n)-amin(n))/interv(n));
        X1 = amin(n) + interv(n)*(ind1-1);
        X2 = amin(n) + interv(n)*(ind2-1);
        dist1 = sqrt((X1-xevo(n))^2);
        dist2 = sqrt((X2-xevo(n))^2);
        if dist1<dist2
            xevor(n) = X1;
        else
            xevor(n) = X2;
        end
    end
[fevor] = predictor(xevor, dmodel); %calls DACE package to obtain predicted function value
next_pts(j,:) = xevor;
min_est(j) = fevor;
end
% next_pts contains several surrogate minimum points, and the surrogate
% poll point (listed last)
disp(strcat('surrogate minimum points are first ', num2str(nsm), ' points in next_pts, surrogate poll point is ', num2str(nsm+1), 'th point listed')) 
if isempty(poll_pts) == 0
    next_pts= [next_pts;poll_xpt]
end
% find MSE max point
% use EA to find maximum MSE point
% minimize -MSE
strfun = 'pred_mse'
xlow = [amin]';
xup = [amax]';
for n = 1:1
    options = es_options_new('TolX',1e-1,'SigmaFacStart',1,'LBound',amin','UBound',amax');
    es = es_new(N, xlow, xup,options);     % generate an ES struct
    es = es_run_mod(es, strfun, 500, dmodel);       % minimize strfun with ES
    result = es_get(es, 'result'); % get result
    xevotmp(:,n) = result.x;               % return current best point
    [fevotmp(n)] = pred_mse(xevotmp(:,n), dmodel); % function pred_mse takes negative of mse to maximize
end
[fevos,I] = sort(fevotmp);
for j = 1:1
    fevo = fevos(j);
    nn = I(j);
    xevo = xevotmp(:,nn);
    for i=1:N
        interv(i) = (amax(i)-amin(i))/(spc(i)-1);
    end
    % find nearest mesh point
    for n=1:N
        ind1 = 1+floor((xevo(n)-amin(n))/interv(n));
        ind2 = 1+ceil((xevo(n)-amin(n))/interv(n));
        X1 = amin(n) + interv(n)*(ind1-1);
        X2 = amin(n) + interv(n)*(ind2-1);
        dist1 = sqrt((X1-xevo(n))^2);
        dist2 = sqrt((X2-xevo(n))^2);
        if dist1<dist2
            xevor(n) = X1;
        else
            xevor(n) = X2;
        end
    end
[fmse(j)] = pred_mse(xevor, dmodel);
max_mse_pt(j,:) = xevor;
end

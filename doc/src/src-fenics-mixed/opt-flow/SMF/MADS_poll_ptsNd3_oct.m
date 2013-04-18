% function to find polling points using MADS, works for rectangular domains, N
% optimization parameters
% finds N+1 poll points which form a positive basis
% no need to switch to 2n near the boundary, since generated poll
% directions are dense
% Alison Marsden  4/27/03
% snap to boundary for the points outside the bounds. 2009.8 wgyang


function poll_pts = MADS_poll_ptsNd3_oct(x,N,delta,spc,amin,amax)

% delta should be 1, 1/4, 1/16, etc.
dg= 1/sqrt(delta);
li = -1/sqrt(delta) + 1;
ui = 1/sqrt(delta) - 1;

rand('state',sum(100*clock))
% fill diag entries with -1 or 1, equal probability
for i=1:N
    if rand < 0.5
        D(i) = dg;
    else
        D(i) = -dg;
    end
end
LT = diag(D);
if N>1
    for i=1:N-1 
        d = diag(LT,-i);
        for j = 1:length(d)
            r = sign(rand-0.5)*(ui)*rand;
            d(j) = round(r);
%             if (3*rand) < 1
%                 d(j) = -1;
%             elseif 1<(3*rand)&(3*rand)<2
%                 d(j) = 0;
%             else
%                 d(j) = 1;
%             end
        end
        LT = LT + diag(d,-i);
    end
end
% randomly permute the rows and columns of LT
rp = randperm(N);
for i = 1:N
    LT1(:,i) = LT(:,rp(i));
end
rp = randperm(N);
for i = 1:N
    LT2(i,:) = LT1(rp(i),:);
end

% complete to an n+1 positive basis
for i=1:N
    interv(i) = (amax(i)-amin(i))/(spc(i)-1);
end
B = [LT2; -sum(LT2)]

poll_pts_pre = [];
%poll_pts1 = [];
%snap to the boundary
for i=1:N+1
    xp = x + interv.*B(i,:); 
    for j=1:N
        if xp(j)>amax(j)
            xp(j)=amax(j);
        end
        if xp(j)<amin(j)
            xp(j)=amin(j);
        end
    end    
        poll_pts_pre = [poll_pts_pre ; xp];
end

num_poll=size(poll_pts_pre,1);

%eliminate repeated points
poll_pts=[];
for i=1:num_poll-1
    repeat=0;
    for j=i+1:num_poll
       d=poll_pts_pre(i,:)-poll_pts_pre(j,:)
       if max(abs(d))<=1e-10    
       repeat=1;
           break
       end
    end
    
    if repeat==0
        poll_pts=[poll_pts ; poll_pts_pre(i,:)];
    end
end

poll_pts=[poll_pts;poll_pts_pre(num_poll,:)];





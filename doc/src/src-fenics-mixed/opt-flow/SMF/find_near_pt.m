% function find_near_pt to find nearest mesh point to a given point

function mesh_pt = find_near_pt(xpt, N, spc, amin, amax)

for i=1:N
        interv(i) = (amax(i)-amin(i))/(spc(i)-1);
    end
for n=1:N
       ind1 = 1+floor((xpt(n)-amin(n))/interv(n));
       ind2 = 1+ceil((xpt(n)-amin(n))/interv(n));
       X1 = amin(n) + interv(n)*(ind1-1);
       X2 = amin(n) + interv(n)*(ind2-1);
       dist1 = sqrt((X1-xpt(n))^2);
       dist2 = sqrt((X2-xpt(n))^2);
       if dist1<dist2
           xpt(n) = X1;
       else
           xpt(n) = X2;
       end
end

mesh_pt = xpt;

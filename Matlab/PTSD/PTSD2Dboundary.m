function [pd, udx, udy] = PTSD2Dboundary(pd, udx, udy, PMLdepth, xiXn, xiXp, xiYn, xiYp)

if xiXn <= 1
    pd(:,PMLdepth+1) = pd(:,PMLdepth+1) .* xiXn;
else
    udx(:,PMLdepth+1) = udx(:,PMLdepth+1) .* (1/xiXn);
end

if xiXp <= 1
    pd(:,end-PMLdepth) = pd(:,end-PMLdepth) .* xiXp;
else
    udx(:,end-PMLdepth) = udx(:,end-PMLdepth) .* (1/xiXp);
end

if xiYn <= 1
    pd(PMLdepth+1,:) = pd(PMLdepth+1,:) .* xiYn;
else
    udy(PMLdepth+1,:) = udy(PMLdepth+1,:) .* (1/xiYn);
end

if xiYp <= 1
    pd(end-PMLdepth,:) = pd(end-PMLdepth,:) .* xiYp;
else
    udy(end-PMLdepth,:) = udy(end-PMLdepth,:) .* (1/xiYp);
end

end

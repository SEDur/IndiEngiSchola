function [pd, udx, udy, udz] = PTSD3Dboundary(pd, udx, udy, udz,...
    PMLdepth, xiXn, xiXp, xiYn, xiYp, xiZn, xiZp)

if xiXn <= 1
    pd(PMLdepth+1,:,:) = pd(PMLdepth+1,:,:) .* xiXn;
else
    udy(PMLdepth+1,:,:) = udy(PMLdepth+1,:,:) .* (1/xiXn);
end

if xiXp <= 1
    pd(end-PMLdepth,:,:) = pd(end-PMLdepth,:,:) .* xiXp;
else
    udy(end-PMLdepth,:,:) = udy(end-PMLdepth,:,:) .* (1/xiXp);
end

if xiYn <= 1
    pd(:,PMLdepth+1,:) = pd(:,PMLdepth+1,:) .* xiYn;
else
    udx(:,PMLdepth+1,:) = udx(:,PMLdepth+1,:) .* (1/xiYn);
end

if xiYp <= 1
    pd(:,end-PMLdepth,:) = pd(:,end-PMLdepth,:) .* xiYp;
else
    udx(:,end-PMLdepth,:) = udx(:,end-PMLdepth,:) .* (1/xiYp);
end

if xiZn <= 1
    pd(:,:,PMLdepth+1) = pd(:,:,PMLdepth+1) .* xiZn;
else
    udz(:,:,PMLdepth+1) = udz(:,:,PMLdepth+1) .* (1/xiZn);
end

if xiZp <= 1
    pd(:,:,end-PMLdepth) = pd(:,:,end-PMLdepth) .* xiZp;
else
    udz(:,:,end-PMLdepth) = udz(:,:,end-PMLdepth) .* (1/xiZp);
end
end

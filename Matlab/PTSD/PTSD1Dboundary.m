function [pd, ud] = PTSD2DboundaryP(pd, ud, PMLdepth, xiL, xiR)
if xiL <= 1
    pd(PMLdepth+1) = pd(PMLdepth+1) * xiL;
else
    ud(PMLdepth+1) = ud(PMLdepth+1) * (1/xiL);
end

if xiR <= 1
    pd(end-PMLdepth) = pd(end-PMLdepth) * xiR;
else
    ud(end-PMLdepth) = ud(end-PMLdepth) * (1/xiR);
end

end

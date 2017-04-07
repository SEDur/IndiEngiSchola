function pd = PTSD3Dsrc(pd, src, srcloc)

pd(srcloc,srcloc,srcloc) = pd(srcloc,srcloc,srcloc) -src;
% pd(srcloc+1,srcloc) = pd(srcloc+1,srcloc) -src;
% pd(srcloc,srcloc+1) = pd(srcloc,srcloc+1) -src;
% pd(srcloc-1,srcloc) = pd(srcloc-1,srcloc) -src;
% pd(srcloc,srcloc-1) = pd(srcloc,srcloc-1) -src;
% pd(srcloc+1,srcloc+1) = pd(srcloc+1,srcloc+1) -src;
% pd(srcloc-1,srcloc-1) = pd(srcloc-1,srcloc-1) -src;
end
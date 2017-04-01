function pd = PTSD2Dsrc(pd, src, srcloc)

pd(srcloc,srcloc) = pd(srcloc,srcloc) -src;
pd(srcloc+1,srcloc) = pd(srcloc+1,srcloc) -src;
pd(srcloc,srcloc+1) = pd(srcloc,srcloc+1) -src;
pd(srcloc-1,srcloc) = pd(srcloc-1,srcloc) -src;
pd(srcloc,srcloc-1) = pd(srcloc,srcloc-1) -src;
pd(srcloc+1,srcloc+1) = pd(srcloc+1,srcloc+1) -src;
pd(srcloc-1,srcloc-1) = pd(srcloc-1,srcloc-1) -src;
end
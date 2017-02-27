function pd = PTSD2Dsrc(pd, src, srcloc, srcgain)

pd(srcloc,srcloc) = pd(srcloc,srcloc) + (src * srcgain);

end
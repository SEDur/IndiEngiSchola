function pd = PTSD2Dsrc(pd, src, srcloc)

pd(srcloc,srcloc) = pd(srcloc,srcloc) + src;

end
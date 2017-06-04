function pd = PTSD3Dsrc(pd, src, srcloc)

% pd(srcloc,srcloc,srcloc) = pd(srcloc,srcloc,srcloc) -src;
pd(srcloc(1)+1,srcloc(2),srcloc(3)+1) = pd(srcloc(1)+1,srcloc(2),srcloc(3)+1) +src;
pd(srcloc(1),srcloc(2)+1,srcloc(3)-1) = pd(srcloc(1),srcloc(2)+1,srcloc(3)-1) +src;
pd(srcloc(1)-1,srcloc(2),srcloc(3)+1) = pd(srcloc(1)-1,srcloc(2),srcloc(3)+1) +src;
pd(srcloc(1),srcloc(2)-1,srcloc(3)-1) = pd(srcloc(1),srcloc(2)-1,srcloc(3)-1) +src;
pd(srcloc(1)+1,srcloc(2)+1,srcloc(3)) = pd(srcloc(1)+1,srcloc(2)+1,srcloc(3)) +src;
pd(srcloc(1)-1,srcloc(2)-1,srcloc(3)) = pd(srcloc(1)-1,srcloc(2)-1,srcloc(3)) +src;
pd(srcloc(1),srcloc(2),srcloc(3)) = pd(srcloc(1),srcloc(2),srcloc(3))+src;

end
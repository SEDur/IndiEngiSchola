function p = FDTD3Dsources(p,sourcelocations ,sourcesamples , sourcetypes)

p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3))...
    = p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3))...
    - sourcesamples(1,1);

end
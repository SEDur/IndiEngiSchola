function [p, ux, uy] = SFDTD2Dfun(p, pCx, pCy, pCt, ux, uy, uCx, uCy, uCt, Rx, Ry, ZL, ZR, ZT, ZB, idx)

if(size(idx(idx > 0),1) > 10)
    %  mat(col, row)
    for i = 1 : size(ux,1)-1
        for i1 = 2 : size(ux,2)-2
            if(idx(i, i1) > 0)
                ux(i,i1) = ux(i,i1) - uCx*(p(i,i1)-p(i,i1-1));
            end
        end
    end
    
    for i = 2 : size(uy,1)-2
        for i1 = 1 : size(uy,2)-1
            if(idx(i, i1) > 0)
                uy(i,i1) = uy(i,i1) - uCx*(p(i,i1)-p(i-1,i1));
            end
        end
    end
else
    % update the non-boundary condition nodes for velocity
    ux(:, 2:end-1) = ux(:, 2:end-1) - uCx*(p(:, 2:end) - p(:, 1:end-1));
    uy(2:end-1, :) = uy(2:end-1, :) - uCy*(p(2:end, :) - p(1:end-1, :));
end

% update the velocity at the right wall
ux(:, end) = ((Rx - ZR)/(Rx + ZR))*ux(:, end) ...
    + (2/(Rx + ZR))*p(:, end);

%update the velocity at the left wall
ux(:, 1) = ((Rx - ZL)/(Rx + ZL))*ux(:, 1) - (2/(Rx + ZL))*p(:, 1);

%update the velocity at the top wall
uy(end, :) = ((Ry - ZT)/(Ry + ZT))*uy(end, :) ...
    + (2/(Ry + ZT))*p(end, :);

%update the velocity at the bottom wall
uy(1, :) = ((Ry - ZB)/(Ry + ZB))*uy(1, :) - (2/(Ry + ZB))*p(1, :);

% update the pressure at all nodes
%     p(idx > 0) = p(idx > 0) - pCx*(ux(idx3 > 0) - ux(idx > 0))...
%         - pCy*(uy(idx2 > 0) - uy(idx > 0));
if(size(idx(idx > 0),1) > 10)
    for i = 1 : size(p,1)-1
        for i1 = 1 : size(p,2)-1
            if(idx(i, i1) > 0)
                p(i,i1) = p(i,i1) - pCx*(ux(i, i1 + 1) - ux(i, i1))...
                    - pCy*(uy(i+1, i1) - uy(i, i1));
            end
        end
    end
else
    p = p - pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
        - pCy*(uy(2:end, :) - uy(1:end-1, :));
end
end
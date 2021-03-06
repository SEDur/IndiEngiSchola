%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD2Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pd, ud] = PSTD2Dfun2(pd, ud, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst, N, Rx, Ry,...
     ZR, ZL, ZF, ZB, pconst, uconst)
    phat = zeros(N,N);
    uhat = zeros(N,N);
    pdiffhat = zeros(N,N);
    udiffhat = zeros(N,N);
    temp = zeros(N, N);
    %% Function solves using the PSTD method for a pressure vector,...
    %  velocity vector and differentiation impulse response in 1 dimension
    %  and returns the solved pressure and velocity vectors
    phat = fft(pd,N,2);
    for i = 1:size(pd, 1)
        temp(i,:) = phat(i,:) .* diffmatrix;
    end
    pdiffhat = ifft(temp,N,2);
%     surf(real(pdiffhat));
    phat = fft(pdiffhat,N,1);
    for i = 1:size(pd, 2)
        temp(:,i) = phat(:,i) .* diffmatrix';
    end
    pdiffhat = ifft(temp,N,1);
%     surf(real(pdiffhat));
    ud = -ud - uconst .* pdiffhat;
    
        % update the velocity at the right wall
    ud(:, end) = ((Rx - ZR)/(Rx + ZR))*ud(:, end) ...
        + (2/(Rx + ZR))*pd(:, end);

    %update the velocity at the left wall
    ud(:, 1) = ((Rx - ZL)/(Rx + ZL))*ud(:, 1)...
        - (2/(Rx + ZL))*pd(:, 1);

    %update the velocity at the top wall
    ud(end, :) = ((Ry - ZF)/(Ry + ZF))*ud(end, :) ...
        + (2/(Ry + ZF))*pd(end, :);

    %update the velocity at the bottom wall
    ud(1, :) = ((Ry - ZB)/(Ry + ZB))*ud(1, :)...
        - (2/(Ry + ZB))*pd(1, :);
    
%     surf(real(ud));
    uhat = fft(ud,N,2);
    for i = 1:size(ud, 1)
        temp(i,:) = uhat(i,:) .* diffmatrix;
    end
    udiffhat = ifft(temp,N,2);
%     surf(real(udiffhat));
    uhat = fft(udiffhat,N,1);
    for i = 1:size(ud, 2)
        temp(:,i) = uhat(:,i) .* diffmatrix';
    end
    udiffhat = ifft(temp,N,1);
%     surf(real(udiffhat));
    pd = -pd - pconst .* udiffhat;
%     surf(real(pd));

end
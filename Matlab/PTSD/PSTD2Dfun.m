%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD1Dfun.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pd, ud] = PSTD2Dfun(pd, ud, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst, N)
    phat = zeros(N,N);
    uhat = zeros(N,N);
    pdiffhat = zeros(N,N);
    udiffhat = zeros(N,N);
    temp = zeros(N, N);
    %% Function solves using the PSTD method for a pressure vector,...
    %  velocity vector and differentiation impulse response in 1 dimension
    %  and returns the solved pressure and velocity vectors
    phat = fft(pd,[],2);
    for i = 1:size(pd, 1)
        temp(i,:) = phat(i,:) .* diffmatrix;
    end
    pdiffhat = ifft(temp,[],2);
    phat = fft(pdiffhat,[],1);
    for i = 1:size(pd, 2)
        temp(:,i) = phat(:,i) .* diffmatrix';
    end
    pdiffhat = ifft(temp,[],1);
    ud = ud .* PMLdiff - PMLalphau .* (pdiffhat./PMLconst); 
    uhat = fft(ud,[],2);
    for i = 1:size(ud, 1)
        temp(i,:) = uhat(i,:) .* diffmatrix;
    end
    udiffhat = ifft(temp,[],2);
    uhat = fft(udiffhat,[],1);
    for i = 1:size(ud, 2)
        temp(:,i) = uhat(:,i) .* diffmatrix';
    end
    udiffhat = ifft(temp,[],1);
    pd = pd .* PMLdiff - PMLalphap .* (udiffhat./PMLconst);

end
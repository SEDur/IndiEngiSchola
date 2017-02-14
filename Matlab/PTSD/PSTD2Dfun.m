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
    for i = 1:size(pd, 1)
        phat(i,:) = fft(pd(i,:));
        temp(i,:) = phat(i,:) .* diffmatrix;
        pdiffhat(i,:) = ifft(temp(i,:));
    end
    for i = 1:size(pd, 2)
        phat(i,:) = fft(pdiffhat(:,1)');
        temp(i,:) = phat(i,:) .* diffmatrix;
        pdiffhat(:,i) = ifft(temp(i,:)');
    end
    ud = ud .* PMLdiff - PMLalphau .* (pdiffhat./PMLconst); 
    for i = 1:size(ud, 1)
        uhat(i,:) = fft(ud(i,:));
        temp(i,:) = uhat(i,:) .* diffmatrix;
        udiffhat(i,:) = ifft(temp(i,:));
    end
    for i = 1:size(ud, 2)
        uhat(i,:) = fft(udiffhat(:,1)');
        temp(i,:) = uhat(i,:) .* diffmatrix;
        udiffhat(:,i) = ifft(temp(i,:)');
    end
    pd = pd .* PMLdiff - PMLalphap .* (udiffhat./PMLconst);

end
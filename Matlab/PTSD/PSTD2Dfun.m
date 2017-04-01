%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTSD2Dfun.m
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
    
%% Velocity in the X dim
%     phat = fft(pd,N,2);
%     for i = 1:size(pd, 1)
%         temp(i,:) = phat(i,:) .* diffmatrix;
%     end
%     pdiffhat = ifft(temp,N,2);
% 
% %% Velocity in the Y Dim
%     phat = fft(pd,N,1);
%     for i = 1:size(pd, 2)
%         temp(:,i) = phat(:,i) .* diffmatrix';
%     end
%     pdiffhat = (pdiffhat + ifft(temp,N,1)) ./2 ;

%% Velocity in 2d
    phat = fft2(pd);
    temp = phat .* diffmatrix;
    pdiffhat = ifft2(temp);

%% Total Velocity
    ud = ud .* PMLdiff - PMLalphau .* (pdiffhat./PMLconst);
    
%% Pressure in 2d
    uhat = fft2(ud);
    temp = uhat .* diffmatrix;
    udiffhat = ifft2(temp);
%% Pressure in the X dim
%     uhat = fft(ud,N,2);
%     for i = 1:size(ud, 1)
%         temp(i,:) = uhat(i,:) .* diffmatrix;
%     end
%     udiffhat = ifft(temp,N,2);
% 
% %% Pressure in the Y Dim
%     uhat = fft(ud,N,1);
%     for i = 1:size(ud, 2)
%         temp(:,i) = uhat(:,i) .* diffmatrix';
%     end
%     udiffhat = (udiffhat + ifft(temp,N,1)) ./2;

%% Total Pressure
    pd = pd .* PMLdiff - PMLalphap .* (udiffhat./PMLconst);


end
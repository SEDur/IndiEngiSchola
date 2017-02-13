% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PTSD1Dfun.m
% % Created by S Durbridge as part of work on a masters dissertation
% % Copywrite S Durbridge 2017
% %
% % Any copies of this function distributed by the autor are done so
% % without any form of warranty, and should not be reproduced without
% % permission
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function[pd, ud] = PSTD1Dfun(pd, ud, diffmatrix,...
%      PMLdiff, PMLalphau, PMLalphap, PMLconst)
%     %% Function solves using the PSTD method for a pressure vector,...
%     %  velocity vector and differentiation impulse response in 1 dimension
%     %  and returns the solved pressure and velocity vectors
%     
%     phat = fft(pd);
%     temp = phat .* diffmatrix;
%     pdiffhat = ifft(temp);
% %     for i2 = 1 : length(pdiffhat)
% %         if i2 < PMLdepth
% %            alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3); 
% %         elseif i2 > N - PMLdepth
% %             alpha = (1/3) * (i2 - ((N-PMLdepth)/PMLdepth)^3);
% %         else
% %             alpha = 0;
% %         end
% %         ud(i2) = ud(i2) * ((1-alpha)/(1+alpha))-uconst * (1/(1+alpha))*(pdiffhat(i2)/(3.142*N));
% %     end
% 
%     ud = ud .* PMLdiff - PMLalphau .* (pdiffhat./PMLconst);
%     
%     uhat = fft(ud);
%     temp = uhat .* diffmatrix;
%     udiffhat = ifft(temp);
%     
% %     for i2 = 1 : length(udiffhat)
% %         if i2 < PMLdepth
% %            alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3); 
% %         elseif i2 > N - PMLdepth
% %             alpha = (1/3) * (i2 - ((N-PMLdepth)/PMLdepth)^3);
% %         else
% %             alpha = 0;
% %         end
% %         pd(i2) = pd(i2) * ((1-alpha)/(1+alpha))-pconst * (1/(1+alpha))*(udiffhat(i2)/(3.142*N));
% %     end
%     pd = pd .* PMLdiff - PMLalphap .* (udiffhat./PMLconst);
% 
% end
function w=hannrectwindow(A,B,C,D,fs)%#codegen
% function w=hannrectwindow(A,B,C,D,[fs])
%
% returns a window where the middle is rectangular and the ends are half
% hann of arbitrary size. A B C D are the indexes of the starts and ends of
% the windows. In principle the rising hann is of length B-A+1 and the falling hann is of
% length D-C+1. The output is always of length of D-A+1.
%
% INPUTS: 
%   A:  index of the pre-hann value of 0(typically A=1)
%   B:  index of the pre-hann value of 1 
%   C:  index of the post-hann value of 1
%   D:  index of the post-hann value of 0 
%   fs: [OPTIONAL] sampling rate. If declared then other args are seconds rather than indexes   
%
% OUTPUTS: w the desired window of length=D-A+1
%
% EXAMPLE:
%   hannrectwindow(1,4,6,12) gives : a 12 point window
%       a (4-1)+1=4 point rising half hann 
%       a (6-4)-1=1 point rect window 
%       a (12-6)+1=7 point falling half hann
%
%   hannrectwindow(1,1,1,12) gives : 
%       a falling half hann of length 12
% 
%
%
% S.Thomson
% B&W Group Ltd, Dec 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=double(A);
B=double(B);
C=double(C);
D=double(D);


if nargin == 5
    % then 1st 4 args are times. convert to samples.
    A =round(A*fs) +1; %offset by one since matlab is not zero bound
    C =round(C*fs) +1;
    B =round(B*fs) +1;
    D =round(D*fs) +1;    
end

if A<=0 || A-B>0 || B-C>0 || C-D>0
    error ('indexes must be positive and in ascending order')
end

if B-A ~=0
    %wpre =hann((B-A+0)*2,'periodic')  ;wpre=wpre(1:B-A+1); %pre window
    wpre =myhann((B-A+0)*2)  ;wpre=wpre(1:B-A+1); %pre window
else
    wpre=[1];
end

if D-C ~=0
    %wpost=hann((D-C+0)*2,'periodic');wpost=flipud(wpost(1:D-C+1)); %post window
    wpost=myhann((D-C+0)*2);wpost=flipud(wpost(1:D-C+1)); %post window
else
    wpost=[1];
end

if C-B ~=0
    wrect=ones(C-B-1,1);
else
    wrect=[];
    %also knock off one of the ones in the pre/post windows
    wpost=wpost(2:end);
end

w=[wpre; wrect; wpost];


function wp=myhann(L)
%return a periodic hann made by hand to support code generation

if coder.target('MATLAB')
    %it's matalb use this
    wp=hann(L,'periodic');
else
    %it's generated code use this
    if L~=1
        Np=L-1+1;
        n=[0:Np-1]';
        wp= 0.5*(1-cos(2*pi*((n)/Np)));
    else
        wp=1;
    end
end

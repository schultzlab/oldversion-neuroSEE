%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Non-Negative Derivative Spike Estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT:
%       R: ratiometric time series whose size is [cells x samples]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spikes = neuroSEE_nnd(R)

% Normalise the Ca traces first
nR = bsxfun(@rdivide, bsxfun(@minus, R, min(R,[],2)),(max(R,[],2)-min(R,[],2)));
nC = size(nR,1);
mAlength = 5;
thr = 0.022; %0.009
spikes = zeros(size(nR));
for ri=1:nC
    nRs = conv(nR(ri,:),ones(mAlength,1)/mAlength,'same');
    temp = diff(nRs);
    temp = [temp,temp(end)];
    temp(temp<thr) = 0;
    spikes(ri,:) = temp;
end

end
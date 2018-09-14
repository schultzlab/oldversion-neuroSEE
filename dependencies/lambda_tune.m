%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambda, phi] = lambda_tune(phi_0, video1, video2, radius,...
                               options, corrIm,...
                               meanIm)
                           
% Settings for tuning lambda
options.corrIm = corrIm;
options.meanIm = meanIm;
options.maxIt = 30;


% Do segmentation
finito      = 0;
lambda(1)   = options.lambda;
counter     = 2;

test_masks = [];
test_masks = squeeze(phi_0(:,:,randi([1,size(phi_0,3)],1,20))); %pick 20 random ROIs from phi_0 
ii = 1;

while ~finito
    [~,~,~,~,~,phi]=segment_cells_tune(test_masks, video1, video2, radius, options);
    
    skew_phi = [];

     
    for k = 1:size(phi,3)
        [x,y] = find(phi(:,:,k)<0);
        temp_phi = [];
        for i=1:length(x)
            for j=1:length(y)
                temp_phi(i) = phi(x(i),y(j),k);
            end
        end
        skew_phi(k) = skewness(temp_phi);
    end
    
    mean_skew_phi = mean(skew_phi);
    
    if mean_skew_phi < 0 % variance is set to minimum of 1
        if length(lambda) == 1
            new_lambda = 0.5*lambda;
        elseif lambda(end) < lambda(end-1)
            new_lambda = 0.5*lambda(end);  
        end
    else 
        finito = 1;
        if new_lambda~=options.lambda %added this ~SEIG
            new_lambda = new_lambda;
        else  
            %new_lambda = lambda; %added this line
            new_lambda = options.lambda; %modified this ~SEIG
        end
    end
    lambda(counter)      = new_lambda;  
    counter              = counter + 1;
    options.lambda       = new_lambda;
    %clearvars A B C D E;
    disp(['Lambda Tuning: Iteration ', num2str(ii), ', current value of lambda is ', num2str(options.lambda)]);
    ii = ii+1;
end
lambda = options.lambda;
disp(['Final lambda is ', num2str(lambda)]);
end
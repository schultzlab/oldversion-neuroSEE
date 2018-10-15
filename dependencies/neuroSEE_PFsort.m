function [placeMap_sorted, alpha_index] = neuroSEE_PFsort(placeMap)

if size(placeMap,2)>360
    disp('Error: cannot perform circular statistics.');
else
    %i = 1:size(placeMap,2);
    i = 0:10:359;
    sin_i = sind(i);
    cos_i = cosd(i);
    alpha = zeros(size(placeMap,1),1); %mean direction
    cv = zeros(size(placeMap,1),1); %circular variance
    for k = 1:size(placeMap,1)
        alpha(k) = atan2d(sum(placeMap(k,:).*sin_i),sum(placeMap(k,:).*cos_i));
        if alpha(k)<0
            alpha(k) = alpha(k)+360;
        end
        cv(k) = 1 - (sqrt(sum(placeMap(k,:).*sin_i)^2 + sum(placeMap(k,:).*cos_i)^2))/360; %circular variance = 1 - ResultantLength
    end
    [~, alpha_index] = sortrows(alpha);
    placeMap_plot = placeMap(alpha_index,:);
    placeMap_sorted=placeMap_plot(:,any(placeMap_plot));
end

end
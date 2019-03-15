% active = getActiveSamples(tracking,thrV)
function active = getActiveSamples(tracking,thrV)
    vx = [0; diff(tracking(:,1)) ];
    vy = [0; diff(tracking(:,2)) ];
    active = sqrt(vx.^2+vy.^2) > thrV;
end
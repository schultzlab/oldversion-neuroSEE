function active = getActiveSamples(tracking,thrV)
    vx = diff(tracking(:,1));
    vy = diff(tracking(:,2));
    active = sqrt(vx.^2+vy.^2)>thrV;
end
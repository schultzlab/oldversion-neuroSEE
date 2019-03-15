function neuroSEE_trackProcess(track_dir,track_names,path_name,sessiondate,mouse);

j=0;
for i=1:length(track_names)
    j = j+1;
    try
        load([track_dir,track_names{i}],'X','Y','phi');
        if mod(size(X,1),2)==0
            x1 = X(2:size(X,1)/2+1,1);
            x2 = X(size(X,1)/2+1:end,1);
            y1 = Y(2:size(Y,1)/2+1,1);
            y2 = Y(size(Y,1)/2+1:end,1);
            xy = [x1 y1];
            phi_1 = phi(2:size(X,1)/2+1,1);
            phi_2 = phi(size(X,1)/2+1:end,1);
            phi = phi_1;
            save([path_name,sessiondate,'_',mouse,'_X',num2str(j),'_track'],'xy','phi', '-v7.3');
            j = j+1;
            xy = [x2 y2];
            phi = phi_2;
            save([path_name,sessiondate,'_',mouse,'_X',num2str(j),'_track'],'xy','phi', '-v7.3');
        else
            x1 = X(2:size(X,1)/2+1,1);
            x2 = X(size(X,1)/2+2:end,1);
            y1 = Y(2:size(Y,1)/2+1,1);
            y2 = Y(size(Y,1)/2+2:end,1);
            xy = [x1 y1];
            phi_1 = phi(2:size(X,1)/2+1,1);
            phi_2 = phi(size(X,1)/2+2:end,1);
            phi = phi_1;
            save([path_name,sessiondate,'_',mouse,'_X',num2str(j),'_track'],'xy','phi', '-v7.3');
            j = j+1;
            xy = [x2 y2];
            phi = phi_2;
            save([path_name,sessiondate,'_',mouse,'_X',num2str(j),'_track'],'xy','phi', '-v7.3');
        end
    catch
        disp(['ERROR: Failed processing data for <', track_names{i}(end-22:end), '>.']);
    end
end
end
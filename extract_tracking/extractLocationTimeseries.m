% [err, track_fname] = extractLocationTimeseries( session, path_name, user_path, track_dir, track_names, logfh )
%
% Extract tracking location time series from csv file & convert to both Cartesian
% & radial formats. Saves a data structure with the following fields:
%    data.tracktime 
%    data.spatialx 
%    data.spatialy   
%    data.spatialphi 

%
% Inputs:
%  session     - date, time, mouse name info in a string
%  path_name   - directory to save tracking info to (on the farm)
%  user_path   - backup user directory if saving to the server fails
%  track_dir   - directory containing the tracking files
%  track_names - cell array of tracking file names (not full path names)
%  logfh       - open file handles that allows writing
% Outputs:
%  err         - error status --> true if successfully saved output tracking file
%  track_fname - name of output tracking file
function [err, track_fname, data] = extractLocationTimeseries( session, path_name, user_path, track_dir, track_names, logfh )
   err = true; track_fname = []; data = [];
   [err, track_names] = neuroSEE_trackProcess( track_dir, track_names, session, path_name, user_path, logfh );
   if err, return; end

   track_time  = [];
   track_xy    = [];
   track_phi   = [];
   track_found = false;
   if ischar( track_names )
      track_names = { track_names };
   end
   for i = 1:length( track_names )
      track_fname = track_names{i};
      if ~exist( track_fname, 'file' )
         str = sprintf( '\tTracking file not found (%s)\n', mc_fname );
         if logfh>0, fprintf( logfh, str ); end
         cprintf('Errors', str);

      else
        track_found = true;
        load( track_fname, 'time', 'xy', 'phi' );
        track_time  = cat(1, track_time, time);
        track_xy    = cat(1, track_xy, xy);
        track_phi   = cat(1, track_phi, phi);
      end
   end
   if ~track_found
      str = sprintf( '\tNo tracking files found, exiting ... \n' );
      if logfh>0, fprintf( logfh, str ); end
      cprintf('Errors', str);
      err = true;
      return
   end

   % Edit the structure names below !!!!
   data.tracktime  = track_time;
   data.spatialx   = track_xy(:,1);
   data.spatialy   = track_xy(:,2);
   data.spatialphi = track_phi;

   % make sure each neuron has a time vector
%    if ~isfield( data, 'imgtime' ) 
%       data.imgtime = 1/scanfreq * (1:length( data.neuron{1}.spikes) );
%    end

   try
      ts_fname = fullfile( path_name, [session '_track.mat'] );
      save( ts_fname, 'data', '-v7.3' ); 
      str = sprintf( '\tTracking masks saved in %s\n', ts_fname );
      if logfh>0, fprintf( logfh, str ); end
      cprintf('Keywords', str);
   catch ME
      % ran out of space so it failed - try saving to my local dir
      if ~exist( user_path, 'dir' )
         mkdir( user_path );
      end
      ts_fname = fullfile( user_path, [session '_track.mat'] );
      save( ts_fname, 'data', '-v7.3' );
      str = sprintf( 'Failed saving tracking file on farm, saving locally (%s)\n', ts_fname );
     if logfh>0, fprintf( logfh, str ); end
      cprintf( 'Errors', str );
   end


end
function groovy_slice(glob_ps, sub_ps)
% slice timing metabatch file
  
% load spm defaults: No need for custom defaults unless you want to
% change the prefix.
spm_get_defaults;

% Ref slice hard coded to 1
ref_slice = 1;

first_flag = 1;
for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  s_filter = ['^' glob_ps.st_prefix this_sub.raw_filter '$'];

  for ss = 1:length(this_sub.sesses) % and session 

	sess = this_sub.sesses(ss);

    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    P = spm_select('List', dirn, s_filter);
    P = [repmat([dirn filesep],size(P,1),1) P];
    
    
    % For 4dnii files
    switch glob_ps.epi_format
        case '4dnii' 
         vol = spm_vol(P);
         filename = vol(1).fname;
		nslices = vol(1).dim(3);
	 sess.acq_order  = [1:2:nslices 2:2:nslices];
        
         for b = 1:length(vol);
             file_list{b}=sprintf('%s, %d',filename,b);
         end;
        
         P = strvcat(file_list);
    end;

    
   % get information from first file, for first session
    if first_flag % first session for each subject
      first_img = deblank(P(1,:));
      V = spm_vol(first_img);
      
      % Sets slice time information
      % value 1 is time to acquire one slice
      % value 2 is time between beginning of last slice
      % and beginning of first slice of next volume
      sl_times = [sess.slice_time ...
		  sess.slice_time + ...
		  (this_sub.TR-sess.slice_time*V.dim(3))];
      firstf = 0;
    end

    % do slice timing correction
    spm_slice_timing(P,sess.acq_order,1,sl_times);

  end
end

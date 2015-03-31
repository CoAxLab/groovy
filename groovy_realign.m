function groovy_realign(glob_ps, sub_ps)

% Get defaults
edefs = glob_ps.realign_opts.estimate;
wdefs = glob_ps.realign_opts.write;
% Flags to pass to routine to calculate realignment parameters
% (spm_realign)
reaFlags = struct(...
    'Quality', edefs.quality,...  % estimation quality
    'fwhm', edefs.fwhm,...        % smooth before calculation
    'sep', edefs.sep,...          % default separation (mm)
    'rtm', edefs.rtm,...          % whether to realign to mean 
    'interp', edefs.interp,...    % interpolation method
    'PW',edefs.weight{1});      % weighting image

% Flags to pass to routine to create resliced images
% (spm_reslice)
resFlags = struct(...
    'interp', wdefs.interp,...    % trilinear interpolation
    'wrap', wdefs.wrap,...        % wrapping info (ignore...)
    'mask', wdefs.mask,...        % masking (see spm_reslice)
    'which',wdefs.which, ...      % whether to write reslice time series
    'mean',1);           % do write mean image

clear imgs; 
% dirnames,  
% get the subdirectories in the main directory

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  r_filter = ['^' glob_ps.realign_prefix this_sub.raw_filter '$'];
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    [P Pdir] = spm_select('List', dirn, r_filter);
    imgs(ss) = {[repmat([dirn filesep],size(P,1),1) P]};
    % For 4dnii files
    switch glob_ps.epi_format
        case '4dnii'
         vol = spm_vol(fullfile(dirn,P));

         try isfield(vol(1),'name');
             filename = vol(1).name;
         catch
             filename = vol(1).fname;
         end;

         for b = 1:length(vol);
             file_list{b}=fullfile(dirn,sprintf('%s, %d',filename,b));
         end;
         
         img(ss) = {strvcat(file_list)};
    end;
    
  end
  
  % Run the realignment
  spm_realign(imgs, reaFlags);
  
  % Run the reslicing
  spm_reslice(imgs, resFlags);
end












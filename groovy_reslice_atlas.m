function groovy_reslice_atlas(glob_ps, sub_ps)

% Get defaults
wdefs = glob_ps.realign_opts.write;

% Flags to pass to routine to create resliced images
% (spm_reslice)
resFlags = struct(...
    'interp', wdefs.interp,...    % trilinear interpolation
    'wrap', wdefs.wrap,...        % wrapping info (ignore...)
    'mask', wdefs.mask,...        % masking (see spm_reslice)
    'which',1, ...      % 1 means leave atlas (first image) the same
    'mean',0);           % do write mean image

clear imgs; 
% dirnames,  
% get the subdirectories in the main directory

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
   r_filter = ['^' glob_ps.smooth_prefix this_sub.raw_filter '$'];
	% war -> rwar prefixes
	
%	clear sess_imgs;
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    [P Pdir] = spm_select('List', dirn, r_filter);
	
	% Add atlas
	atlasfile = fullfile('/','data','templates-atlases','AAL625.nii');
	imgs(1,1) = {atlasfile};
    imgs(2,1) = {[repmat([dirn filesep],size(P,1),1) P]};
		% Builds NIFTI filename and saves it to second row
  
  % Run the reslicing
  spm_reslice(imgs, resFlags); % imgs is really just one image
	end
end

function groovy_norm_calc(glob_ps, sub_ps)
% Normalization metabatch; calculates parameters and writes structural only 

% Get the defaults
defs = glob_ps.normalise_opts;

% Turn of template weighting
defs.estimate.weight = '';

for s = 1:length(sub_ps) % for each subject 
  my_sub = sub_ps(s);

	for ss = 1:length(my_sub.sesses) % for each session
		sess = my_sub.sesses(ss);
  
	  % Get normalization source image
	  subj_sess(ss).P = sess.norm_source;
	  %dirn = fullfile(glob_ps.fdata_root, my_sub.dir);
	  %pfile = spm_select('List', dirn, ['^' my_sub.raw_filter]);
	  %subj(s).P = [dirn filesep 'mean' pfile];
	  
	  % Make the default normalization parameters file name
	  subj_sess(ss).matname = ...
		[spm_str_manip(subj_sess(ss).P,'sd') '_sn.mat'];
	  
	  % Set the object mask for subject
	  subj_sess(ss).objmask = my_sub.obj_mask;
	
	  % set orientation (.mat file) just in case
	  if ~isempty(subj_sess(ss).objmask) 
	    spm_get_space(subj_sess(ss).objmask, ...
			  spm_get_space(subj_sess(ss).P));
	  end
	  
	  % Get the images we are going to reslice
	  % Because we are going reslice later 
	  % We don't reslice anything except the image to be normalized
	  subj_sess(ss).PP = subj_sess(ss).P;      
	  % call the SPM normalize function to do the work
	  spm_normalise(glob_ps.template_images, ...
			subj_sess(ss).P, subj_sess(ss).matname,...
			defs.estimate.weight, subj_sess(ss).objmask, ...
			defs.estimate);
	  
	  % Do the reslicing
	  spm_write_sn(subj_sess(ss).PP,subj_sess(ss).matname,defs.write);
	end
end

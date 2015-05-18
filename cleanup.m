function cleanup(glob_ps, sub_ps)

% {global_params.fdata_root}/{subject}/{rs_dir}/{sub_struct.raw_filter}.nii

for s = 1:length(sub_ps) % for each subject
	this_sub = sub_ps(s);
	for ss = 1:length(this_sub.sesses) % for each session 
		dirn = fullfile(glob_ps.fdata_root, ...
				this_sub.dir, this_sub.sesses(ss).dir);
		
		% Select original and preprocessed file
		filter = ['^(rwar)?' this_sub.raw_filter '$'];
		allfiles = spm_select('List', dirn);

		for i = 1:size(allfiles, 1)
			filename = strtrim(allfiles(i,:));
			if isempty(regexp(filename, filter))
				filepath = strtrim([dirn filesep filename]);
				delete(filepath)
				%msg = sprintf('Deleted %s', filepath);
				%disp(msg)
			end
		end
	end
end

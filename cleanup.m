function cleanup(glob_ps, sub_ps)

% {global_params.fdata_root}/{subject}/{rs_dir}/{sub_struct.raw_filter}.nii

for s = 1:length(sub_ps) % for each subject
	this_sub = sub_ps(s);
	for ss = 1:length(this_sub.sesses) % for each session 
		dirn = fullfile(glob_ps.fdata_root, ...
				this_sub.dir, this_sub.sesses(ss).dir);
		
		% Select original and preprocessed file
		if glob_ps.gradunwarp
			n_filters = 3
		else
			n_filters = 2;
		end
		filters = cell(n_filters, 1);
		filters{1} = ['^(rwar)?' this_sub.raw_filter '$'];
		filters{2} = ['^rp_.*\.txt$'];
		if glob_ps.gradunwarp
			filters{3} = ['^wmean.*$']
		end
		allfiles = spm_select('List', dirn);

		for i = 1:size(allfiles, 1)
			filename = strtrim(allfiles(i,:));
			matches_filter = 0;
			for idx=1:n_filters
				if ~isempty(regexp(filename, filters{idx}))
					matches_filter = 1;
				end
			end
			if ~matches_filter
				filepath = strtrim([dirn filesep filename]);
				delete(filepath);
				%msg = sprintf('Deleted %s', filepath);
				%disp(msg)
			end
		end
	end
end

function uncompress(glob_ps, sub_ps)

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
	% Select all gzipped NIFTIs
   r_filter = ['^' glob_ps.realign_prefix this_sub.raw_filter '.gz$'];

  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    [P Pdir] = spm_select('List', dirn, r_filter);
	
keyboard
	gunzip(P)
	
	end
end

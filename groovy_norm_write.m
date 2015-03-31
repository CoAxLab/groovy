function groovy_norm_write(glob_ps, sub_ps)
% metabatch file to write with normalization parameters for SPM8

defs = glob_ps.normalise_opts;

clear imgs;
for s = 1:length(sub_ps) % for each subject 
  my_sub = sub_ps(s);
  subj_dir = fullfile(glob_ps.fdata_root, my_sub.dir);
  
  % Make the default normalization parameters file name
  ns = my_sub.norm_source; 
  matname = [spm_str_manip(ns, 'sd') '_sn.mat'];
  
  if iscell(my_sub.norm_others);
      my_sub.norm_others = char(my_sub.norm_others);
  end;
  
  % Do the reslicing for miscellaneous images
  spm_write_sn(my_sub.norm_others,matname,defs.write);

  % Do the reslicing for the EPIs
  nw_filter = ['^' glob_ps.norm_write_prefix my_sub.raw_filter '$'];
  imgs = '';
  for ss = 1:length(my_sub.sesses)
    dirn = fullfile(subj_dir,my_sub.sesses(ss).dir);
    P = spm_select('List',dirn,nw_filter);
    
    switch glob_ps.epi_format
      case '4dnii';
        vol = spm_vol(fullfile(dirn,P));
        filename = vol(1).fname;
        for b = 1:length(vol);
          file_list{b}=sprintf('%s, %d',filename,b);
        end;
        imgs = strvcat(imgs,strvcat(file_list));
        
      otherwise
        imgs = strvcat(imgs, [repmat([dirn filesep],size(P,1),1) P]);
    end;
    
  end
  spm_write_sn(imgs,matname,defs.write);
end



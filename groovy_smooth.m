function groovy_smooth(glob_ps, sub_ps)

% Get the defaults
sdefs = glob_ps.smooth_opts;

imgs = '';
for sb = 1:length(sub_ps)
  my_sub = sub_ps(sb);
  subjroot = fullfile(glob_ps.fdata_root, my_sub.dir);
  sm_filter = ['^' glob_ps.smooth_prefix my_sub.raw_filter '$'];
  for ss = 1:length(my_sub.sesses)
    dirn = fullfile(subjroot,my_sub.sesses(ss).dir);
    % get files in this directory
    P = spm_select('List',dirn,sm_filter);
    
    % For 4dnii files
    switch glob_ps.epi_format
        case '4dnii' 
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
end
  
% and this is just spm_smooth_ui pasted/edited
s     = sdefs.fwhm;
P     = imgs;
n     = size(P,1);
dtype = sdefs.dtype;

% implement the convolution
%---------------------------------------------------------------------------
for i = 1:n
  Q = deblank(P(i,:));
  [pth,nm,xt] = fileparts(deblank(Q));
  U = fullfile(pth,[sdefs.prefix nm xt]);
  spm_smooth(Q,U,s);
end





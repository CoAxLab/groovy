function groovy_diagnostics(glob_ps, sub_ps)
% slice timing metabatch file
  
% load spm defaults: No need for custom defaults unless you want to
% change the prefix.
spm_get_defaults;

% Ref slice hard coded to 1
ref_slice = 1;

first_flag = 1;
for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  s_filter = ['^' glob_ps.diag_prefix this_sub.raw_filter '$'];
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    P = spm_select('List', dirn, s_filter);
    P = [repmat([dirn filesep],size(P,1),1) P];
    
    
    % For 4dnii files
    switch glob_ps.epi_format
        case '4dnii' 
         vol = spm_vol(P);
         filename = vol(1).fname;
        
         for b = 1:length(vol);
             file_list{b}=sprintf('%s, %d',filename,b);
         end;
        
         P = strvcat(file_list);
    end;

    V = spm_vol(P);
    Y = spm_read_vols(V);
    
    % Get mean & stdev images
    mY = squeeze(mean(Y,4));
    sY = squeeze(std(Y,[],4));
    
    mV = V(1); mV.fname = fullfile(glob_ps.fdata_root, this_sub.dir, ...
				   this_sub.sesses(ss).dir,'meanEPI.nii');
    
    sV = V(1); sV.fname = fullfile(glob_ps.fdata_root, this_sub.dir, ...
				   this_sub.sesses(ss).dir,'stdEPI.nii');

    s2nV = V(1); s2nV.fname = fullfile(glob_ps.fdata_root, this_sub.dir, ...
				   this_sub.sesses(ss).dir,'s2nEPI.nii');
    
    spm_write_vol(mV, mY);
    spm_write_vol(sV, sY);
    spm_write_vol(s2nV, mY./sY);
    
    figure;
    mid_vox = round(size(mY,3)/2);
    subplot(1,3,1);
    imagesc(squeeze(mY(:,:,mid_vox))); colorbar; axis square
    subplot(1,3,2);
    imagesc(squeeze(sY(:,:,mid_vox))); colorbar; axis square
    subplot(1,3,3);
    
    s2nY = mY./sY;
    imagesc(squeeze(s2nY(:,:,mid_vox))); colorbar; axis square;
    
    saveas(gcf,fullfile(glob_ps.fdata_root, ...
			this_sub.dir,this_sub.sesses(ss).dir, ...
			'meanEPI_plots.fig'));

    close;
    
    % Look at head motion
    rp_fname = dir(fullfile(glob_ps.fdata_root, this_sub.dir, ...
			    this_sub.sesses(1).dir,'rp_*.txt'));
    rp_file = fullfile(glob_ps.fdata_root, this_sub.dir, ...
		       this_sub.sesses(1).dir,rp_fname.name);
    mov = load(rp_file);
    
    figure;
    subplot(2,1,1);
    plot(mov(:,1:3));
    legend('X','Y','Z');
    ylabel('Shift (mm)');
    xlabel('Time (vol)')
    title('Translation');
    
    subplot(2,1,2);
    plot(mov(:,4:6));
    legend('Pitch','Roll','Yaw');   
    ylabel('Degrees');
    xlabel('Time (vol)')    
    title('Rotation');
    saveas(gcf,fullfile(glob_ps.fdata_root, this_sub.dir, ...
     			 this_sub.sesses(1).dir, ...
     			 'motion_params.fig'));
    
    
     close;

  end
end

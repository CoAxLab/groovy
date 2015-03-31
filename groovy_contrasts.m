function groovy_contrasts(glob_ps, sub_ps)
% batch file to run contrasts

if ~exist('defaults', 'var'), spm_get_defaults; end

% store path
pwd_orig = pwd;

for s = 1:length(sub_ps) % for each subject 
  this_sub = sub_ps(s);
  
  % get, goto SPM results directory
  ana_dir = fullfile(glob_ps.fdata_root, ...
		     this_sub.dir, ...
		     glob_ps.ana_sdir);
  cd(ana_dir);
  
  % load SPM model; give "SPM" structure
  disp('Loading SPM.mat');
  load('SPM.mat');
  disp('Done');
  
  % Fix swd, just in case
  SPM.swd = ana_dir;
  
  % Where we will start filling in the contrasts
  if ~isfield(SPM,'xCon');
      SPM.xCon = [];
  end;
  
  xcon_1 = length(SPM.xCon)+1;
  
  % now put contrast into SPM structure
  for cn = 1:length(this_sub.contrasts)
      struct = spm_FcUtil('Set',...
          this_sub.contrasts(cn).name,...
          this_sub.contrasts(cn).type,...
          'c', ...
          this_sub.contrasts(cn).con_mat',...
          SPM.xX.xKXs);
      
      if ~size(SPM.xCon,2) 
          SPM.xCon = struct;
      else
          SPM.xCon(end + 1) = struct;
      end;
      
  end;
  
  % Estimate only the contrasts we've added
  spm_contrasts(SPM, xcon_1:length(SPM.xCon));
  
end

% back to initial directory
cd(pwd_orig);

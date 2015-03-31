function groovy_subject_model(glob_ps, sub_ps)
% metabatch file to run SPM2 models for single subjects
global defaults
spm_get_defaults

spm fmri;

% store path
pwd_orig = pwd;
   
% Specify design stuff common to all subjects
%===========================================================================
% global normalization: OPTOINS:'Scaling'|'None'
%----------------------------------------------------------spm_select('List',fildir,'-----------------
gSPM.xGX.iGXcalc       = 'None';

% low frequency confound: high-pass cutoff (secs) [Inf = no filtering]
%---------------------------------------------------------------------------
gSPM.xX.K(1).HParam    = glob_ps.high_pass;

% intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w'
%-----------------------------------------------------------------------
%
%SPM.xVi.form       = 'AR(1) + w';
gSPM.xVi.form       = glob_ps.ar_form;

% basis functions and timing parameters
%---------------------------------------------------------------------------
% OPTIONS:'hrf'
%         'hrf (with time derivative)'
%         'hrf (with time and dispersion derivatives)'
%         'Fourier set'
%         'Fourier set (Hanning)'
%         'Gamma functions'
%         'Finite Impulse Response'
%---------------------------------------------------------------------------
% Fill in the field below with the corresponding string above
gSPM.xBF.name       = glob_ps.event_bf.name;

% length in seconds - not used for hrf 
gSPM.xBF.length     = glob_ps.event_bf.length; 

% order of basis set - not used for hrf
gSPM.xBF.order      = glob_ps.event_bf.order;  

% The next two fields usually don't need changing. 
% number of time bins per scan
gSPM.xBF.T          = 16;      

% first time bin (see slice timing
gSPM.xBF.T0         = 1;       

% Selfish explanatory - OPTIONS: 'scans'|'secs' for	onsets
gSPM.xBF.UNITS      = 'scans';

% value of one means no Volterra pain - OPTIONS: 1|2 = order of
% convolution
gSPM.xBF.Volterra   = 1;                 

for sb = 1:length(sub_ps)
  this_sub = sub_ps(sb);
  
  % Get template SPM structure
  SPM = gSPM;
  
  % Note that the TR must be the same for all runs in a model
  SPM.xY.RT          =  this_sub.TR;    % seconds
  
  % specify filter for filenames
  Filter             = ['^' glob_ps.stats_prefix this_sub.raw_filter '$'];
  
  % get, make, goto SPM results directory
  sub_dir = fullfile(glob_ps.fdata_root,this_sub.dir);
  ana_sdir = fullfile(sub_dir, glob_ps.ana_sdir);
  if ~(exist(ana_sdir))
    mkdir(sub_dir,glob_ps.ana_sdir);
  end
  cd(ana_sdir);
  
  PP=[];
  
  for ss = 1:length(this_sub.sesses)
    % Information for this session
    this_ss = this_sub.sesses(ss);
    
    % directory containing scans
    fildir = fullfile(sub_dir, this_ss.dir);
    
    % file selection
    P = spm_select('List',fildir,Filter);
    
    % For 4dnii files
    switch glob_ps.epi_format
        case '4dnii'
         vol = spm_vol(fullfile(fildir,P));
         filename = vol(1).fname;
         
         clear file_list
         for b = 1:length(vol);
             file_list{b}=sprintf('%s, %d',filename,b);
         end;
         P = strvcat(file_list)
    end;

    n_scans = size(P,1);
    SPM.nscan(ss) = n_scans;
    PP = strvcat(PP, P);
    
    % Condition stuff - onsets, durations, types.
    for cno = 1:length(this_ss.cond_names)
      ons = this_ss.ons{cno};
      dur = this_ss.dur{cno};

      try
        pmod = this_ss.pmod(cno);
	if isempty(pmod.name);
	  pmod = struct('name','none');
	end;

      catch
        pmod = struct('name','none');
      end;

      early_ons = ons < 0;
      if any(early_ons)
       warning('Some onsets are less than 0');
       disp(ons(early_ons))
     end
     
     late_ons = ons > n_scans;
     if any(late_ons), 
       warning('Some onsets are after the end of scanning');
       disp(ons(late_ons))
     end

     ons = ons(~(early_ons | late_ons));
     if isempty(ons)
       error(...
         sprintf('Session %d, %d scans; no onsets are left',...
          ss, n_scans));
     end
     dur = dur(~(early_ons | late_ons));

     SPM.Sess(ss).U(cno) = struct(...
       'ons', ons, ...
       'dur', dur,...
       'name',{this_ss.cond_names(cno)},...
	     'P', pmod); % Parametric modulation
   end

    % Movement stuff
    movefil = spm_select('List',fildir,'^rp.*txt$');
    %movefil = spm_get('Files', fildir, ['rp_*.txt']);
    [m_c m_names] = movement_regressors(fullfile(fildir,movefil), glob_ps.movement_params);
    
    if ~isempty(m_c)
      % Fix any goofy realignment params
      m_c = m_c(1:n_scans, :);
    end
    
    % other covariates and names go in first
    covs = [this_ss.covs m_c];
    cov_names = [this_ss.cov_names m_names];
    
    % Put covariates into model
    SPM.Sess(ss).C.C    = covs;     % [n x c double] covariates
    SPM.Sess(ss).C.name = cov_names; % [1 x c cell]   names
    
  end % session loop
  
  % set files
  SPM.xY.P           = PP;
  
  % Configure design matrix
  SPM = spm_fmri_spm_ui(SPM);
  
  % Estimate parameters
  spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
  
  % Estimate the model using the appropriate algorithm
  if ~isfield(glob_ps,'algorithm')
      SPM = spm_spm(SPM);
  else
      switch glob_ps.algorithm
          case 'spm'
              SPM = spm_spm(SPM);
          case 'rwls'
              SPM = spm_rwls_spm(SPM);
          otherwise
              error('Unknown Model Estimation Type in groovy_subject_model');
      end;
  end;
  
end

% back to initial directory
cd(pwd_orig);

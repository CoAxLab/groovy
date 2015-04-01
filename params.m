function [global_params, subject_params] = params(subjects)

% sets up various things specific to this analysis
% You'll be using your own version of this function, so its name will
% change.  In the comments below, the function name will be
% "this_function"
% 
% FORMAT [global_params, subject_params] = this_function(subjects)
%   
% Input
% subjects       - subject directory name identifying subject OR
%                  cell array of subject directory names OR
%                  vector of subject numbers
%                  If empty, all subjects' data returned
% 
% Outputs
% global_params  - parameters shared across all subjects' analyses
% subject_params - struct array of parameters per subject
%
% The fields in global_params contain all the stuff that is common for all
% subjects - like the session names.
%
% The subject_params struct array contains one struct per subject, with one
% field ('dir') giving the subject directory name - e.g subject_params(1) =
% struct('dir', '00AH', 'sone_things', 1, 'other_things', 2) and so on.
% This way we dont necessarily have to keep track of things like subject
% _numbers_ we can reference the subjects details with a string - e.g. to
% get the TR for a particular subject:
% 
% [gp sp] = this_function('00AH');
% TR = sp.TR
% 
% OR
%  
% [gp sp] = this_function;
% sub_no = strmatch('00AH', {sp(:).dir}, 'exact');
% TR = sp(sub_no).TR
% 
% First written by MB, 30 March 2005
% Modified for SPM8 by T. Verstynen, 18 Feb 2010
% Modified for SPM12 by M Miller, 30 March 2015

% Where the subjects' data directories are stored
% Defaults to current working directory
global_params.fdata_root = pwd;

% Can uncomment to set the fdata_root by computer
%switch computer
%	case 'MACI64'
%     		global_params.fdata_root = '/Users/tmoleswo/Data/WIN/';
%	otherwise
%		global_params.fdata_root = pwd;
%end;

all_subjects = {'hcp'};
    
if nargin < 1
  subjects = [];
end

if isempty(subjects)
  subjects = all_subjects;
elseif isnumeric(subjects)
  subjects = all_subjects(subjects);
elseif ischar(subjects)
  subjects = cellstr(subjects); 
end




% Where the directory tree containing parameters is stored 
% This could be where you have stored your batch files, some already
% calculated normalization parameters, reference functions and so on.
% parameter_root could be the% same as fdata_root above.
global_params.parameter_root = fullfile(global_params.fdata_root, ...
    'groovy_batch');

% Whether to store structural stuff on parameter or fdata root
my_anat_root = global_params.parameter_root;

% Session directories, in fact same for each subject
my_sesses = {...
    'BOLD_resting_PMU',...
    };


nsubs = length(subjects);
nsesses = length(my_sesses);
nslices = 72;
ntrs = 210;

% prefixes for processed images (prepended to raw image filter below, in
% the scripts)
global_params.st_prefix = 'r';  
global_params.realign_prefix = ''; 
global_params.norm_write_prefix = 'ar';
global_params.diag_prefix = 'r';
global_params.smooth_prefix = 'war';
global_params.stats_prefix = 'ar';
global_params.restingstate_prefix = 'swar';

% condition stuff - in fact same for each subject and session 
cond_names = {}; 

con_mat = {};
con_names = {};

con_type = {};
	    
% covariate names - in fact same for each subject and session
% For none, need empty cell array "{}"
cov_names = {};

% Length of trials in TRs
trial_dur = 1;


% ***************************************************
% Slice Time Correction
% ***************************************************
% Because things happen when this changes across subjects, we need to
% have this defined above so as to push the values into each the subject
% structure.

% Slice Time Duration: 
global_params.slicetime = 0.010;

% For Interleaved Acquisition
acq_order = [1:2:nslices-1 2:2:nslices];

% ***************************************************
% Physio info
% ***************************************************

target_physio = '';  % options 'ECG', 'RESP','PPOLOW','PPOHIGH', 'NONE';

% Now set subject by subject stuff
for sb = 1:nsubs
  sub_str = subjects{sb};
  sub_dir_f = fullfile(global_params.fdata_root, sub_str);
  
  % Filter for raw image names.
  rs_dir = 'connport';
  sub_struct.raw_filter = ['100307_fnca_BOLD_REST2_LR.nii'];
  
  % image to normalize for this subject
  sub_struct.norm_source = fullfile(sub_dir_f,rs_dir,...
      sprintf('mean%s',sub_struct.raw_filter));
  
  sub_struct.n_trs = ntrs;
  
  
  % Other images (in space of structural) to write normalized
  % (epis resliced in anothe write-normalized pass)
  sub_struct.norm_others = fullfile(sub_dir_f,rs_dir,...
      [global_params.stats_prefix sub_struct.raw_filter]);
  
  %sub_struct.norm_others = fullfile(sub_dir_f,'rest1', ...    
  % 				    sub_str, ...
  % 				    'anatomical.img');
  
  % object mask - usually empty
  sub_struct.obj_mask = '';
  
  % TR for each subject.  Sometimes it's different for each subject
  % but in this case it's the same
  sub_struct.TR = 0.720;

  % Contrast information for each subject
  sub_struct.con_mat = con_mat;
  sub_struct.con_names = con_names;
  
  for ss = 1:nsesses
    % Session structure, within subject structure
    ss_struct = [];
        
    % Fill session structure
    % here the directory names are all the same
    ss_struct.dir = [rs_dir]; 
    
     % Condition file    
    ss_struct.ons = {};
    ss_struct.dur = {};
    
    % Contrast information for each subject
%     sub_struct.contrasts(c).name = con_names{c};
%     sub_struct.contrasts(c).type = upper(con_type{c});
%     sub_struct.contrasts(c).con_mat = con_mat{c};
%   
    ss_struct.cond_names = cond_names;        
    ss_struct.covs = [];
    ss_struct.cov_names = cov_names;
   
    switch target_physio
     case 'NONE'
      physio_parms = '';
     case 'ECG+RESP'
      physio_params = {'ECG','RESP'};
     otherwise
      physio_params = {target_physio}; 
    end;
    
    if ~isempty(target_physio)
      phlem_file = spm_select('List',fullfile(sub_dir_f,ss_struct.dir),...
			      '^PhysioPhLEM.*mat');
      
      if isempty(phlem_file)
	fprintf(sprintf('Unable to find PhLEM output file in dir:\n \t %s \n RECONSTRUCTING \n',ss_struct.dir))
	
	physio_f = dir(fullfile(sub_dir_f,ss_struct.dir,'PhysioLog*.*'));
	ecg_file = fullfile(sub_dir_f,ss_struct.dir,physio_f(1).name);
	resp_file = fullfile(sub_dir_f,ss_struct.dir,physio_f(3).name);
	pulse_file = fullfile(sub_dir_f,ss_struct.dir,physio_f(2).name);
	
	
	PhLEM_setup(fullfile(sub_dir_f,ss_struct.dir,'PhysioPhLEM.mat'),...
		    'ecgfile',ecg_file,'respfile',resp_file,...
		    'ppofile',pulse_file);
	% PhLEM_setup(ecg_file,resp_file,pulse_file);
	phlem_file = spm_select('List',fullfile(sub_dir_f,ss_struct.dir),...
				'^PhysioPhLEM.*mat');
      end;
      
      tmp = load(fullfile(sub_dir_f,ss_struct.dir,phlem_file));        
      PhLEM = tmp.PhLEM;
      
      for ph = 1:length(physio_params);    
	[C names] = make_physio_regressors(physio_params{ph},PhLEM,...
					   'TR',sub_struct.TR);
    
	ss_struct.covs = [ss_struct.covs C];
	ss_struct.cov_names = {ss_struct.cov_names{:} names{:}};
      end;
      
    end;  
    
    % Put into subject structure
    sub_struct.sesses(ss) = ss_struct;
  end

  % Directory and filter for coreg target image
  sub_struct.coreg_target_dir = fullfile(sub_dir_f, sub_struct.sesses(1).dir);
  sub_struct.coreg_target_filter = ['mean' sub_struct.raw_filter];
  
  % And to-be-coregistered image
  sub_struct.coreg_object_dir = fullfile(my_anat_root, sub_str);
  sub_struct.coreg_object_filter = 'bs*mpflash.img';
  
  % And images to take along with coregistration
  % dirs field here has to be cell array to cope with many dirs, for
  % sessions
  sub_struct.coreg_other_dirs = {};
  sub_struct.coreg_other_filter = 'bs*mpflash.img';
    
  % Set the subdirectory for this subject.  This is required 
  sub_struct.dir = sub_str;
  
  % Put the slice time in for this subject (Edit here
  % if an error occurs and it is different per subject
  sub_struct.slice_time = global_params.slicetime;
  sub_struct.acq_order = acq_order;
  
  % Set into returned structure
  subject_params(sb) = sub_struct;
  
end

% Setup the Global pointers to use in key preprocessing routines
global defaults
spm_get_defaults;


% ***************************************************
% Motion Correction and Realignment 
% ***************************************************
% Estimation Parameters: see spm_realign for more
%
%  QUALITY: Quality versus speed trade-off.  Highest quality
%           (1) gives most precise results.  Default = 0.9
% 
%  SEP: the default separation (mm) to sample the images.
%
%  RTM: register to mean. Default = 1
% 
%  PW: a filename of a weighting image (reciprocal of
%      standard deviation)
% 
%  INTERP: B-spline degree used for interpolation. Def = 2;
%

% Set the fields.
global_params.realign_opts.estimate = defaults.realign.estimate;

% Modify to realign to the first image instead of mean
global_params.realign_opts.estimate.rtm = 0;

% Modify to doing a simple trilinear interpolation (1st order)
global_params.realign_opts.estimate.interp = 2;

% Don't use a wegithing image
global_params.realign_opts.estimate.weight = {[]};


% Writing Parameters
%
% PREFIX: prefix for resliced images. Defaults to 'r'.
% 
% MASK: mask output images (1 for yes, 0 for no)
%
% MEAN: write mean image (1 = yes, 0 = no)
%
% INTERP: the B-spline interpolation method. 
% 
% WHICH: what to output from the realignment (def = 1 2);
%       0   - don't create any resliced images.
%             Useful if you only want a mean resliced image.
%       1   - don't reslice the first image.
%       2   - reslice all the images.

% Set the fields.
global_params.realign_opts.write = defaults.realign.write;

% Change the reslicing default to all but first.
global_params.realign_opts.write.which = [2];


% ***************************************************
% Smoothing
% ***************************************************
% 
% PREFIX: prefix for resliced images. Defaults to 's'.
%
% FWHM: [sx, sy, sz] Gaussian filter width in mm
%
% DTYPE: datatype [[default: 0 == same datatype as P]

% Set the defaults:
global_params.smooth_opts = defaults.smooth;

% Smoothing for normalized images
global_params.smooth_opts.fwhm = repmat(4,1,3);

% Smoothing for normalized images
global_params.smooth_opts.dtype = 0;

% ***************************************************
% Normalization: Recommend manual parameter estimation
% ***************************************************
%
% Estimation Parameters: Template weighting is turned off here
%
% SMOSRC: smoothing of source image (FWHM of Gaussian in mm).
%          Defaults to 8.
% SMOREF: smoothing of template image (defaults to 0).
% REGTYPE: regularisation type for affine registration
%           See spm_affreg.m (default = 'mni').
% CUTOFF: Cutoff of the DCT bases.  Lower values mean more
%          basis functions are used (default = 30mm).
% NITS: number of nonlinear iterations (default=16).
% REG:  amount of regularisation (default=0.1)

% Set the defaults
global_params.normalise_opts = defaults.old.normalise;

% Change the source image smoothing to 4mm FWHM
global_params.normalise_opts.estimate.smorc = 4;

%
% Writing Parameters: 
%
%  PREFIX: prefix for resliced images. Defaults to 'w'
%
%  WRAP: wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%
%  VOX: voxel sizes (3 element vector - in mm)
%       Non-finite values mean use template vox.
%       (default - 2x2x2mm)
%
%  BB: bounding box (2x3 matrix - in mm)
%     Non-finite values mean use template bb.
%
%  PRESERVE: either 0 or 1.  A value of 1 will "modulate"
%            the spatially normalised images so that total
%            units are preserved, rather than just
%            concentrations.
%
%  INTERP:  interpolation method (0-7; Default 1-- Trilinear)

% Change the voxel sizes to 1mm
global_params.normalise_opts.vox = [1 1 1];


% List of templates for normalization
% Here, a smoothed grey matter segmentation of MNI brain 
% Have to determine which version of SPM we're using
% in order to determine where the template images are
global_params.template_images = ... 
	fullfile(spm('Dir'), 'toolbox', 'OldNorm', 'EPI.nii');


% ***************************************************
% GLM Parameters: 
% ***************************************************

% subdirectory name for analysis
global_params.ana_sdir = 'RS_GLM';

% Event basis function
% OPTIONS:'hrf'
%         'hrf (with time derivative)'
%         'hrf (with time and dispersion derivatives)'
%         'Fourier set'
%         'Fourier set (Hanning)'
%         'Gamma functions'
%         'Finite Impulse Response'
global_params.event_bf.name = 'hrf';
global_params.event_bf.length = 20;  % this isn't used for an hrf model
global_params.event_bf.order = 1;    % this isn't used for an hrf model

% High pass filter in seconds (Inf= no filtering)
global_params.high_pass = Inf;

% intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w'
global_params.ar_form = 'AR(1)+w';

% What movement parameters to inlude in the model
% Cell array of one or more of 
% 'moves'   - user Nx^ movement parameters as regressors
% 'mm1'     - use movement parameters lagged by 1
% 'moves_2' - movement parameters squared
% 'mm1_2'   - mm1 squared
% (empty cell array = no movement parameter correction)
global_params.movement_params = {'moves'};

% Which estimation routine to use.
% Cell array of one of
% 'spm'   - use the standard SPM model (Default if this isn't present)
% 'rwls'  - use the reweighted least-squares toolbox
global_params.algorithm = 'rwls';

% What is the general data format of the EPIs?
% Cell array of one of
% 'analyze' - 3-D analyze images (.hdr/.img)
% '3dnii'   - 3-D nifti files
% '4dnii'   - 4-D nifti files 
global_params.epi_format = '4dnii';


% ***************************************************
% MELODIC Parameters: 
% ***************************************************

% The bandpass filters for the butterworth filter [High Low]
global_params.bp_filter = [31.25 1.667];

% The new prefix for the filtered file name
global_params.bp_prefix = 'bp';

% Voxel threshold for mask image.  Only selects intensities above this
% voxel value
global_params.mask_thresh = 600;

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
% The full path to a sample unprocessed, uncompressed
% .nii file is:
% {global_params.fdata_root}/{subject}/{rs_dir}/{sub_struct.raw_filter}.nii
% subjects can be provided or chosen automatically through filters
% Be sure to check:
%  * rs_dir
%  * sub_struct.raw_filter

global_params.fdata_root = pwd;

% Can uncomment to set the fdata_root by computer
%switch computer
%	case 'MACI64'
%     		global_params.fdata_root = '/Users/tmoleswo/Data/WIN/';
%	otherwise
%		global_params.fdata_root = pwd;
%end;

% If no subject names are specified, defaults to all the directories within
% global_params.fdata_root, filtered by subj_filter, regular expression
global_params.subj_filter = 'sub\d+';

listing = dir(global_params.fdata_root);
dirs = {listing([listing.isdir]==1).name}; % Get all directories
for idx = length(dirs):-1:1 % Filter directories by subj_filter
	if isempty(regexp(dirs{idx}, 'sub\d+'))
		dirs(idx) = [];
	end
end
    
if nargin < 1
  subjects = [];
end

if isempty(subjects)
  subjects = dirs;
elseif isnumeric(subjects)
  subjects = all_subjects(subjects);
elseif ischar(subjects)
  subjects = cellstr(subjects); 
end

nsubs = length(subjects);

% Scan-specific information
% Be sure to review and edit:
%  * nslices
%  * ntrs
%  * global_params.slicetime (TR / nslices)
%  * sub_struct.TR (see below)
nslices = 30;
% Check number of slices using spm_vol

ntrs = 210;

% ***************************************************
% Slice Time Correction
% ***************************************************
% Because things happen when this changes across subjects, we need to
% have this defined above so as to push the values into each the subject
% structure.

% Slice Time Duration: 
%global_params.slicetime = 0.051;

% For Interleaved Acquisition
acq_order = [1:2:nslices-1 2:2:nslices];

% Where the directory tree containing parameters is stored 
% This could be where you have stored your batch files, some already
% calculated normalization parameters, reference functions and so on.
% parameter_root could be the% same as fdata_root above.
global_params.parameter_root = fullfile(global_params.fdata_root, ...
    'groovy_batch');

% Whether to store structural stuff on parameter or fdata root
my_anat_root = global_params.parameter_root;

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
% Physio info
% ***************************************************

target_physio = '';  % options 'ECG', 'RESP','PPOLOW','PPOHIGH', 'NONE';

% Now set subject by subject stuff
for sb = 1:nsubs
  sub_str = subjects{sb};
  sub_dir_f = fullfile(global_params.fdata_root, sub_str);

  % rs_dir is where within the subject directory data is stored
  % Separate multiple directories with a '/'
	% my_sesses will fill with all the subdirectories in subject dirs
  rs_dir = 'BOLD';

  % Session directories: all directories within subject folders.
	% If that isn't the case, simply change rs_dir to something
	% other than my_sesses and set my_sesses = ''

	%listing = dir([sub_dir_f])
	% Comment above line and uncomment below if my_sesses is beneath rs_dir
  listing = dir([sub_dir_f filesep rs_dir]);
  my_sesses = {listing([listing.isdir]==1).name};
  my_sesses(1:2) = []; % remove '.' and '..' directories

   nsesses = length(my_sesses);

  % The filter for files within rs_dir
  % Make sure the file name includes this pattern
  sub_struct.raw_filter = ['bold.nii'];

  
  % TR for each subject.
  sub_struct.TR = 1.54;
  sub_struct.n_trs = ntrs;
  
  %sub_struct.norm_others = fullfile(sub_dir_f,'rest1', ...    
  % 				    sub_str, ...
  % 				    'anatomical.img');

  % object mask - usually empty
  sub_struct.obj_mask = '';
  
  % Contrast information for each subject
  sub_struct.con_mat = con_mat;
  sub_struct.con_names = con_names;

  % Directory and filter for coreg target image
  	sub_struct.coreg_target_dir = sub_dir_f;
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
  sub_struct.dir = [sub_str filesep rs_dir];

	if isfield(sub_struct, 'sesses')
		sub_struct= rmfield(sub_struct, 'sesses'); % clear sesses struct
	end

 	% If have sessions 
    % Session structure, within subject structure
	if nsesses
		for ss = 1:nsesses
		  ss_struct = [];
		      
		  % Fill session structure
		  % here the directory names are all the same
		  %ss_struct.dir = [my_sesses(ss)]; 
		  ss_struct.dir = my_sesses{ss}; 
			filepath = fullfile(sub_dir_f, rs_dir, ss_struct.dir);

		  % Find the file that fits the filter
			%if class(rs_dir) ~= 'char'
			  %pfile = spm_select('List', fullfile(sub_dir_f, ss_struct.dir), ...
				  %['^' sub_struct.raw_filter]);
			%else
				pfile = spm_select(...
							'List', filepath,['^' sub_struct.raw_filter '$']);
				compressed = spm_select('List', filepath,...
								 ['^' sub_struct.raw_filter '.gz']);

			if size(pfile, 1) > 1
				disp(sprintf('%i NIFTI file found in a session. Only 1 allowed'))
			end
			
			% Uncompress NIFTIs if compressed
			for idx = 1:size(compressed, 1)
				filename = compressed(idx, :);
				filepath = fullfile(sub_dir_f, rs_dir, ss_struct.dir, ...
										filename);
				[pathstr name ext] = fileparts(filename);
				if strcmp(ext, '.gz')
					already_unzipped = 0;
					for pidx = 1:size(pfile, 1) % Check if already unzipped
						if ~strcmp([pfile(1,pidx) '.gz'], compressed(1,:))
							already_unzipped = 1;
						end
					end
					if ~already_unzipped
						gunzip(filepath);
						dirpath = fullfile(sub_dir_f, rs_dir, ss_struct.dir);
						pfile = spm_select(...
							'List', dirpath,['^' sub_struct.raw_filter '$']);
					end
				end
			end

		  % image to normalize for this subject
		  ss_struct.norm_source = fullfile(sub_dir_f,rs_dir, ss_struct.dir,...
		    sprintf('mean%s',pfile));  

			% Other images (in space of structural) to write normalized
			% (epis resliced in another write-normalized pass)
			ss_struct.norm_others = fullfile(sub_dir_f, rs_dir, ss_struct.dir,...
			    [global_params.stats_prefix pfile]);

		   % Condition file    
		  ss_struct.ons = {};
		  ss_struct.dur = {};
		  
		  % Contrast information for each subject
%		    sub_struct.contrasts(c).name = con_names{c};
%		    sub_struct.contrasts(c).type = upper(con_type{c});
%		    sub_struct.contrasts(c).con_mat = con_mat{c};
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

		      % Calculate slice information with SPM from the file
			filepath = fullfile(sub_dir_f, rs_dir, ss_struct.dir, pfile);
			[pathstr name ext] = fileparts(filepath);
			%if strcmp(ext, '.gz')
			%	gunzip(filepath);
			%	filepath = [pathstr name];
		      volume_test = spm_vol(filepath);
		      nslices = volume_test(1).dim(3);
		
		       ss_struct.slice_time = sub_struct.TR./nslices;

		       % for interleaved acquisition
			% rewritten in groovy_slice.m from realigned image
			ss_struct.acq_order  = [1:2:nslices 2:2:nslices];
		  
		  % Put into subject structure
		  sub_struct.sesses(ss) = ss_struct;
		end

	else % Make default session variables for when have no sessions
		ss_struct.dir = '';

		  % Find the file that fits the filter
		  pfile = spm_select('List', fullfile(sub_dir_f, rs_dir), ...
		      ['^' sub_struct.raw_filter '$']);
		compressed = spm_select('List', filepath,...
							 ['^' sub_struct.raw_filter '.gz']);

		% Calculate slice information with SPM from the file
			filepath = fullfile(sub_dir_f, rs_dir, ss_struct.dir, pfile);
			[pathstr name ext] = fileparts(filepath)

			% Uncompress NIFTIs if compressed
			for idx = 1:size(compressed, 1)
				filename = compressed(idx, :);
				filepath = fullfile(sub_dir_f, rs_dir, ss_struct.dir, ...
										filename);
				[pathstr name ext] = fileparts(filename);
				if strcmp(ext, '.gz')
					already_unzipped = 0;
					for pidx = 1:size(pfile, 1) % Check if already unzipped
						if ~strcmp([pfile(1,pidx) '.gz'], compressed(1,:))
							already_unzipped = 1;
						end
					end
					if ~already_unzipped
						gunzip(filepath);
						pfile = spm_select(...
							'List', filepath,['^' sub_struct.raw_filter '$']);
					end
				end
			end
		volume_test = spm_vol(filepath);
		nslices = volume_test(1).dim(3);
		
		ss_struct.slice_time = sub_struct.TR./nslices;

		% for interleaved acquisition
		% rewritten in groovy_slice.m from realigned image
		ss_struct.acq_order  = [1:2:nslices 2:2:nslices];

		  % image to normalize for this subject
		  ss_struct.norm_source = fullfile(sub_dir_f,ss_struct.dir,...
		    sprintf('mean%s',pfile));  

			% Other images (in space of structural) to write
			% normalized (epis resliced in anothe write-normalized 
			% pass)
			ss_struct.norm_others = fullfile(sub_dir_f, ss_struct.dir,...
			    [global_params.stats_prefix pfile]);

		sub_struct.sesses(1) = ss_struct; % put session structure in
	end
  
  % Put the slice time in for this subject (Edit here
  % if an error occurs and it is different per subject
  %sub_struct.slice_time = global_params.slicetime;
  %sub_struct.acq_order = acq_order;

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

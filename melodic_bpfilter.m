function PP = melodic_bpfilter(glob_ps, sub_ps)
% metabatch file to run SPM2 models for single subjects

% Physical pointer to FSL code
switch computer
    case 'MACI64'
     fsl_str='!/Applications/fsl/bin/fslmaths';
    
    otherwise
     fsl_str = '!/usr/bin/fsl5.0-fslmaths';
end;

% Find BP filter ranges
filter_range = glob_ps.bp_filter;

% Get the output prefix
oprefix = glob_ps.bp_prefix;

% store path
pwd_orig = pwd;

PP = [];

for sb = 1:length(sub_ps)
  this_sub = sub_ps(sb);
    
  % specify filter for filenames
  Filter             = ['^' glob_ps.restingstate_prefix this_sub.raw_filter '$'];
  
  % get, make, goto SPM results directory
  sub_dir = fullfile(glob_ps.fdata_root,this_sub.dir);
    
  for ss = 1:length(this_sub.sesses)
    % Information for this session
    this_ss = this_sub.sesses(ss);
    
    % directory containing scans
    fildir = fullfile(sub_dir, this_ss.dir);
    
    % file selection
    P = spm_select('List',fildir,Filter);
    
    % Input output files
    i_file = fullfile(fildir,P);
    o_file = fullfile(fildir,sprintf('%s%s',oprefix, P));
    
    % Form the command
    eval_str = sprintf('%s %s -bptf %2.3f %2.3f %s -odt float',...
        fsl_str, i_file, filter_range(1), filter_range(2), o_file);
    
    % Evaluate the command
    eval(eval_str);

    % Save the file strings
    PP = strvcat(PP, o_file);
  end
end;
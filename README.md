# groovy
ConnPort preprocessing for resting-state fMRI data.

Groovy is a number of MATLAB scripts using SPM functions to pre-process
resting-state fMRI data.

## Requirements
* MATLAB r2014a (tested)
* SPM12

## Installation
Download and unzip the scripts to a directory, then be sure to add that 
directory to MATLAB's search path. You can test installation by typing 
`which batch_preprocess` in MATLAB and checking the file path.

## Configuration
Groovy takes unprocessed, uncompressed NIFTI EPI files as input.

The **params.m** file has settings to adjust. It's best to copy or move
params.m to the parent directory of subject data. Then you can save
different settings for different scans as separate params.m files.

### Specifying the input data directory
In **params.m**, the complete file path to your data will be:

`{global_params.fdata_root}/{subject}/{rs_dir}/{sub_struct.raw_filter}.nii`

The variable `global_params.fdata_root` will by default be the MATLAB current 
working directory. You can also specify a different path in params.m.

The `all_subjects` variable contains a list of subject IDs. The script will
look for those IDs as separate directories when running. By default,
`all_subjects` is populated with all the directories in the 
`global_params.fdata_root`.

Within those subject ID directories, the script will look into `rs_dir`, which
can be a filepath with multiple directories (just separate with the usual `/`).
By default, this is set to the `my_sesses` variable, which by default contains
all the directories within each subject directory. If you don't have separate
sessions, set `my_sesses = '';` and `rs_dir` to whatever directory inside the
subject directory the data is in. Set `rs_dir = '';` if the images are within
the subject directories themselves. The value of `rs_dir` can specify longer
filepaths, for example, `rs_dir = '<top_level_dir>/<lower_level_dir>'`.

Finally, the script has a regular expression filter, `sub_struct.raw_filter`, 
to find input files within `rs_dir` directories. Construct a regular expression
pattern that will match all unprocessed files you wish to preprocess.

### Specify other scan-specific parameters
Other important parameters to review in **params.m** include
* TR: `sub_struct.TR`
* number of volumes: `ntrs`

The number of slices per volume is calculated per image (using `spm_vol`). The
slice time duration is calculated as TR/(number of slices per volume).

## Running
After configuring **params.m**, simply run `batch_preprocess` in MATLAB to 
preprocess the data. (`batch_preprocess` needs to be in MATLAB's search 
path).

## Output
Preprocessed files are placed in the same directory as the input data. Prefixes
are attached to the filenames to signify that they are the output of a certain
preprocessing step.

The prefixes:
* `r{filename}` slice-time corrected
* `ar{filename}` realigned
* `war{filename} normalized to MNI space

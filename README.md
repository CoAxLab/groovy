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
params.m to the parent directory of subject data. That way you can save
different settings for different scans as separate params.m files.

### Specifying the input data directory
In **params.m**, the complete file path to your data will be:

`{global_params.fdata_root}/{subject}/{rs_dir}/{sub_struct.raw_filter}.nii`

The variable `global_params.fdata_root` will by default be the MATLAB current 
working directory. You can also specify a different path in params.m.

The `all_subjects` variable contains a list of subject IDs. The script will
look for those IDs as separate directories when running.

Within those subject ID directories, the script will look into `rs_dir`, which
can be a filepath with multiple directories (just separate with the usual `/`).

Finally, the script has a filter, `sub_struct.raw_filter`, to find files within
that `rs_dir` directory. Make your unprocessed file will pass that filter
(having the same name works).

### Specify other scan-specific parameters
Other important parameters to review in **params.m** include
* TR: `sub_struct.TR`
* number of slices per volume: `nslices`
* slice time (TR/(number of slices per volume)): `global_params.slicetime`
* number of volumes: `ntrs`

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

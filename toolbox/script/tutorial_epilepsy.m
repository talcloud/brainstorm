function tutorial_epilepsy(tutorial_dir)
% TUTORIAL_EPILEPSY: Script that reproduces the results of the online tutorial "EEG/Epilepsy".
%
% CORRESPONDING ONLINE TUTORIALS:
%     http://neuroimage.usc.edu/brainstorm/Tutorials/Epilepsy
%     http://neuroimage.usc.edu/brainstorm/Tutorials/EpilepsyScript
%
% INPUTS: 
%     tutorial_dir: Directory where the sample_epilepsy.zip file has been unzipped

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Author: Francois Tadel, 2014


% ===== FILES TO IMPORT =====
% You have to specify the folder in which the tutorial dataset is unzipped
if (nargin == 0) || isempty(tutorial_dir) || ~file_exist(tutorial_dir)
    error('The first argument must be the full path to the tutorial dataset folder.');
end
% Build the path of the files to import
AnatDir   = fullfile(tutorial_dir, 'sample_epilepsy', 'anatomy');
RawFile   = fullfile(tutorial_dir, 'sample_epilepsy', 'data', 'tutorial_eeg.bin');
ElcFile   = fullfile(tutorial_dir, 'sample_epilepsy', 'data', 'tutorial_electrodes.elc');
SpikeFile = fullfile(tutorial_dir, 'sample_epilepsy', 'data', 'tutorial_spikes.txt');
% Check if the folder contains the required files
if ~file_exist(RawFile)
    error(['The folder ' tutorial_dir ' does not contain the folder from the file sample_epilepsy.zip.']);
end

% ===== CREATE PROTOCOL =====
% The protocol name has to be a valid folder name (no spaces, no weird characters...)
ProtocolName = 'TutorialEpilepsy';
% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 1);


% ===== BRAINSTORM-GENERATED SCRIPT =====
% Start a new report
bst_report('Start');
% Subject name
SubjectName = 'sepi01';

% ===== ANATOMY =====
% Process: Import anatomy folder
bst_process('CallProcess', 'process_import_anatomy', [], [], ...
    'subjectname', SubjectName, ...
    'mrifile',     {AnatDir, 'FreeSurfer'}, ...
    'nvertices',   15000, ...
    'nas', [135, 222,  75], ...
    'lpa', [ 57, 118,  68], ...
    'rpa', [204, 119,  76], ...
    'ac',  [131, 145, 110], ...
    'pc',  [130, 119, 111], ...
    'ih',  [128, 134, 170]);

% Process: Generate BEM surfaces
bst_process('CallProcess', 'process_generate_bem', [], [], ...
    'subjectname', SubjectName, ...
    'nscalp',      1922, ...
    'nouter',      1922, ...
    'ninner',      1922, ...
    'thickness',   4);


% ===== LINK CONTINUOUS FILE =====
% Process: Create link to raw file
sFilesRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
    'subjectname',    SubjectName, ...
    'datafile',       {RawFile, 'EEG-DELTAMED'}, ...
    'channelreplace', 1, ...
    'channelalign',   0);

% Process: Set channel file
sFilesRaw = bst_process('CallProcess', 'process_import_channel', sFilesRaw, [], ...
    'channelfile',  {ElcFile, 'XENSOR'}, ...
    'usedefault',   1, ...
    'channelalign', 1);

% Process: Set channels type
sFilesRaw = bst_process('CallProcess', 'process_channel_settype', sFilesRaw, [], ...
    'sensortypes', 'SP1, SP2, RS, PHO, DELR, DELL, QR, QL', ...
    'newtype',     'MISC');

% Process: Project electrodes on scalp
sFilesRaw = bst_process('CallProcess', 'process_channel_project', sFilesRaw, []);

% Process: Events: Import from file
sFilesRaw = bst_process('CallProcess', 'process_evt_import', sFilesRaw, [], ...
    'evtfile', {SpikeFile, 'ARRAY-TIMES'}, ...
    'evtname', 'SPIKE');

% Process: Snapshot: Sensors/MRI registration
bst_process('CallProcess', 'process_snapshot', sFilesRaw, [], ...
    'target',   1, ...  % Sensors/MRI registration
    'modality', 4, ...  % EEG
    'orient',   1, ...  % left
    'comment',  'MEG/MRI Registration');


% ===== EVALUATION 50 Hz =====
% Process: Power spectrum density (Welch)
sFilesPsd = bst_process('CallProcess', 'process_psd', sFilesRaw, [], ...
    'timewindow',  [], ...
    'win_length',  10, ...
    'win_overlap', 50, ...
    'clusters',    {}, ...
    'sensortypes', 'EEG', ...
    'edit', struct(...
         'Comment',         'Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'none', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));

% Process: Snapshot: Frequency spectrum
bst_process('CallProcess', 'process_snapshot', sFilesPsd, [], ...
    'target',   10, ...  % Frequency spectrum
    'modality', 4, ...   % EEG
    'comment',  'Power spectrum density');


% ===== HIGH-PASS FILTER =====
% Process: Import MEG/EEG: Time
sFilesRaw = bst_process('CallProcess', 'process_import_data_time', sFilesRaw, [], ...
    'subjectname', SubjectName, ...
    'condition',   '', ...
    'timewindow',  [], ...
    'split',       0, ...
    'usectfcomp',  1, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% Process: High-pass:0.5Hz
sFilesRaw = bst_process('CallProcess', 'process_bandpass', sFilesRaw, [], ...
    'highpass',    0.5, ...
    'lowpass',     0, ...
    'mirror',      1, ...
    'sensortypes', '', ...
    'overwrite',   1);

% Process: Compute noise covariance
bst_process('CallProcess', 'process_noisecov', sFilesRaw, [], ...
    'baseline', [120, 130], ...
    'dcoffset', 1, ...
    'method',   1, ...  % Full noise covariance matrix
    'copycond', 0, ...
    'copysubj', 0);


% ===== IMPORT EVENTS ===
% Process: Import MEG/EEG: Events
sFilesEpochs = bst_process('CallProcess', 'process_import_data_event', sFilesRaw, [], ...
    'subjectname', SubjectName, ...
    'condition',   '', ...
    'eventname',   'SPIKE', ...
    'timewindow',  [], ...
    'epochtime',   [-0.1, 0.3], ...
    'createcond',  1, ...
    'ignoreshort', 1, ...
    'usectfcomp',  1, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);

% Process: Average: By condition (subject average)
sFilesAvg = bst_process('CallProcess', 'process_average', sFilesEpochs, [], ...
    'avgtype',    3, ...
    'avg_func',   1, ...  % Arithmetic average: mean(x)
    'keepevents', 0);

% Process: Snapshot: Recordings time series
bst_process('CallProcess', 'process_snapshot', sFilesAvg, [], ...
    'target',   5, ...  % Recordings time series
    'modality', 4, ...  % EEG
    'comment',  'Average spike');

% Process: Snapshot: Recordings topography (contact sheet)
bst_process('CallProcess', 'process_snapshot', sFilesAvg, [], ...
    'target',   7, ...  % Recordings topography (contact sheet)
    'modality', 4, ...  % EEG
    'orient',   1, ...  % left
    'contact_time',   [-0.040, 0.110], ...
    'contact_nimage', 16, ...
    'comment',  'Average spike');

% Process: Import MEG/EEG: Time
sFilesBaseline = bst_process('CallProcess', 'process_import_data_time', sFilesRaw, [], ...
    'subjectname', SubjectName, ...
    'condition',   '', ...
    'timewindow',  [120, 130], ...
    'split',       0, ...
    'usectfcomp',  1, ...
    'usessp',      1, ...
    'freq',        [], ...
    'baseline',    []);


% ===== SOURCE MODELING =====
% Process: Compute head model
bst_process('CallProcess', 'process_headmodel', sFilesAvg, [], ...
    'sourcespace', 1, ...
    'eeg',         3, ...  % OpenMEEG BEM
    'openmeeg',    struct(...
         'BemSelect',    [0, 0, 1], ...
         'BemCond',      [1, 0.0125, 1], ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemFiles',     {{}}, ...
         'isAdjoint',    0, ...
         'isAdaptative', 1, ...
         'isSplit',      0, ...
         'SplitLength',  4000));

% Process: Compute sources
sFilesSrc = bst_process('CallProcess', 'process_inverse', [sFilesAvg, sFilesBaseline], [], ...
    'method', 1, ...  % Minimum norm estimates (wMNE)
    'wmne', struct(...
         'NoiseCov',      [], ...
         'InverseMethod', 'wmne', ...
         'ChannelTypes',  {{}}, ...
         'SNR',           3, ...
         'diagnoise',     0, ...
         'SourceOrient',  {{'fixed'}}, ...
         'loose',         0.2, ...
         'depth',         1, ...
         'weightexp',     0.5, ...
         'weightlimit',   10, ...
         'regnoise',      1, ...
         'magreg',        0.1, ...
         'gradreg',       0.1, ...
         'eegreg',        0.1, ...
         'ecogreg',       0.1, ...
         'seegreg',       0.1, ...
         'fMRI',          [], ...
         'fMRIthresh',    [], ...
         'fMRIoff',       0.1, ...
         'pca',           1), ...
    'sensortypes', 'EEG', ...
    'output', 1);  % Kernel only: shared

% Process: Snapshot: Sources (one time)
bst_process('CallProcess', 'process_snapshot', sFilesSrc(2), [], ...
    'target',   8, ...  % Sources (one time)
    'modality', 1, ...  % MEG (All)
    'orient',   3, ...  % top
    'time',     0, ...
    'comment', 'Average spike');


% ===== Z-SCORE =====
% Process: Z-score normalization: [120.000s,130.000s]
sFilesZscore = bst_process('CallProcess', 'process_zscore_ab', sFilesSrc(1), sFilesSrc(2), ...
    'baseline',   [120, 130], ...
    'source_abs', 1);

% Process: Snapshot: Sources (one time)
bst_process('CallProcess', 'process_snapshot', sFilesZscore, [], ...
    'target',   8, ...  % Sources (one time)
    'modality', 1, ...  % MEG (All)
    'orient',   3, ...  % top
    'time',     0, ...
    'comment', 'Average spike');

% Save and display report
ReportFile = bst_report('Save', sFilesZscore);
bst_report('Open', ReportFile);

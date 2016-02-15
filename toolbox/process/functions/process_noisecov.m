function varargout = process_noisecov( varargin )
% PROCESS_NOISECOV: Compute a noise covariance matrix

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
% Authors: Francois Tadel, 2012-2016

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'Compute noise covariance';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 321;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NoiseCovariance';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Option: Baseline
    sProcess.options.baseline.Comment = 'Time window:';
    sProcess.options.baseline.Type    = 'baseline';
    sProcess.options.baseline.Value   = [];
    % Option: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG, SEEG, ECOG';
    % Option: noisecov/ndatacov
    sProcess.options.label0.Comment = '<BR>Matrix to estimate:';
    sProcess.options.label0.Type    = 'label';
    sProcess.options.target.Comment = {'Noise covariance', 'Data covariance'};
    sProcess.options.target.Type    = 'radio';
    sProcess.options.target.Value   = 1;
    % Options: Remove DC offset
    sProcess.options.label1.Comment = 'Remove DC offset:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.dcoffset.Comment = {'Block by block, to avoid effects of slow shifts in data', 'Compute global average and remove it to from all the blocks'};
    sProcess.options.dcoffset.Type    = 'radio';
    sProcess.options.dcoffset.Value   = 1;
    % Option: Full/Diagonal
    sProcess.options.method.Comment = {'Full noise covariance matrix', 'Diagonal matrix (better if: nTime < nChannel*(nChannel+1)/2)', 'No noise modeling (identity matrix)'};
    sProcess.options.method.Type    = 'radio';
    sProcess.options.method.Value   = 1;
    sProcess.options.method.Group   = 'Output';
    % Option: Copy to other conditions
    sProcess.options.copycond.Comment = 'Copy to other conditions';
    sProcess.options.copycond.Type    = 'checkbox';
    sProcess.options.copycond.Value   = 0;
    sProcess.options.copycond.Group   = 'Output';
    % Option: Copy to other subjects
    sProcess.options.copysubj.Comment = 'Copy to other subjects';
    sProcess.options.copysubj.Type    = 'checkbox';
    sProcess.options.copysubj.Value   = 0;
    sProcess.options.copysubj.Group   = 'Output';
    % Option: Replace file
    sProcess.options.replacefile.Comment = {'Replace', 'Merge', 'Keep', 'If file already exists: '};
    sProcess.options.replacefile.Type    = 'radio_line';
    sProcess.options.replacefile.Value   = 1;
    sProcess.options.replacefile.Group   = 'Output';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    isDataCov = (sProcess.options.target.Value == 2);
    % Get default options
    OPTIONS = bst_noisecov();
    % Get options
    isIdentity = 0;
    if isfield(sProcess.options, 'baseline') && isfield(sProcess.options.baseline, 'Value') && iscell(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value{1})
        OPTIONS.Baseline = sProcess.options.baseline.Value{1};
    else
        OPTIONS.Baseline = [];
    end
    if isfield(sProcess.options, 'sensortypes') && ~isempty(sProcess.options.sensortypes)
        OPTIONS.ChannelTypes = strtrim(str_split(sProcess.options.sensortypes.Value, ','));
    else 
        OPTIONS.ChannelTypes = [];
    end
    switch (sProcess.options.dcoffset.Value)
        case 1,  OPTIONS.RemoveDcOffset = 'file';
        case 2,  OPTIONS.RemoveDcOffset = 'all';
    end
    switch (sProcess.options.method.Value)
        case 1,  OPTIONS.NoiseCovMethod = 'full';
        case 2,  OPTIONS.NoiseCovMethod = 'diag';
        case 3,  isIdentity = 1;
    end
    % Copy to other studies
    isCopyCond = sProcess.options.copycond.Value;
    isCopySubj = sProcess.options.copysubj.Value;
    % Replace file?
    if isfield(sProcess.options, 'replacefile') && isfield(sProcess.options.replacefile, 'Value') && ~isempty(sProcess.options.replacefile.Value)
        switch (sProcess.options.replacefile.Value)
            case 1,   OPTIONS.ReplaceFile = 1;
            case 2,   OPTIONS.ReplaceFile = 2;
            case 3,   OPTIONS.ReplaceFile = 0;   
        end
    else
        OPTIONS.ReplaceFile = 1;
    end

    % ===== GET DATA =====
    % Get all the input data files
    iStudies = [sInputs.iStudy];
    iDatas   = [sInputs.iItem];
    % Get channel studies
    [tmp, iChanStudies] = bst_get('ChannelForStudy', iStudies);
    % Keep only once each channel file
    iChanStudies = unique(iChanStudies);
    
    % ===== COMPUTE =====
    % No noise modeling: Use identity matrix
    if isIdentity
        NoiseCovFiles = import_noisecov(iChanStudies, 'Identity', OPTIONS.ReplaceFile, isDataCov);
    % Compute NoiseCov matrix
    else
        NoiseCovFiles = bst_noisecov(iChanStudies, iStudies, iDatas, OPTIONS, isDataCov);
    end
    if isempty(NoiseCovFiles)
        bst_report('Error', sProcess, sInputs, 'Unknown error.');
        return;
    end
        
    % ===== GET OUTPUT STUDY =====
    % Only the input studies
    if ~isCopyCond && ~isCopySubj
        iCopyStudies = [];
    % All the conditions of the selected subjects
    elseif isCopyCond && ~isCopySubj
        iCopyStudies = [];
        AllSubjFile = unique({sInputs.SubjectFile});
        for iSubj = 1:length(AllSubjFile)
            [tmp, iNew] = bst_get('StudyWithSubject', AllSubjFile{iSubj});
            iCopyStudies = [iCopyStudies, iNew];
        end
    % The selected conditions for all the subjects
    elseif ~isCopyCond && isCopySubj
        iCopyStudies = [];
        ProtocolSubjects = bst_get('ProtocolSubjects');
        AllCond = unique({sInputs.Condition});
        AllSubj = {ProtocolSubjects.Subject.Name};
        for iSubj = 1:length(AllSubj)
            for iCond = 1:length(AllCond)
                [tmp, iNew] = bst_get('StudyWithCondition', [AllSubj{iSubj}, '/', AllCond{iCond}]);
                iCopyStudies = [iCopyStudies, iNew];
            end
        end
    % All the studies
    elseif isCopyCond && isCopySubj
        ProtocolStudies = bst_get('ProtocolStudies');
        iCopyStudies = 1:length(ProtocolStudies.Study);
    end
    iCopyStudies = unique(iCopyStudies);
    
    % ===== COPY TO OTHER STUDIES =====
    if ~isempty(iCopyStudies)
        % Get channel studies
        [tmp, iCopyChanStudies] = bst_get('ChannelForStudy', iCopyStudies);
        % Remove studies that are already processed
        iCopyChanStudies = setdiff(unique(iCopyChanStudies), iChanStudies);
        % Copy noise covariance to other subjects/conditions (overwrites)
        if ~isempty(iCopyChanStudies)
            db_set_noisecov(iChanStudies(1), iCopyChanStudies, isDataCov, OPTIONS.ReplaceFile);
        end
    end
    % Return the data files in input
    OutputFiles = {sInputs.FileName};
end




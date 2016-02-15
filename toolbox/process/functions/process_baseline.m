function varargout = process_baseline( varargin )
% PROCES_BASELINE: Remove the baseline average from each channel (for the given time instants).
%
% DESCRIPTION: For each channel:
%   1) Compute the mean m for the baseline
%   2) For all the time samples, subtract m

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
% Authors: Francois Tadel, 2010-2015

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Remove DC offset';
    sProcess.FileTag     = '| bl';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 60;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/Epoching?highlight=%28Remove+DC+offset%29#Import_in_database';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data', 'results', 'matrix'};
    sProcess.OutputTypes = {'raw', 'data', 'results', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Default values for some options
    sProcess.processDim  = 2;    % Process time by time

    % Definition of the options
    sProcess.options.description.Comment = ['For each signal in input:<BR>' ...
                                            '1) Compute the mean <I>m</I> for the baseline<BR>' ...
                                            '2) For each time sample, subtract <I>m</I><BR><BR>'];
    sProcess.options.description.Type    = 'label';
    % === Baseline time window
    sProcess.options.baseline.Comment = 'Baseline:';
    sProcess.options.baseline.Type    = 'baseline';
    sProcess.options.baseline.Value   = [];
    % === Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data', 'raw'};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    % Get baseline
    if isfield(sProcess.options, 'baseline') && isfield(sProcess.options.baseline, 'Value') && iscell(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value{1})
        time = sProcess.options.baseline.Value{1};
    else
        time = [];
    end
    % Comment: seconds or miliseconds
    if isempty(time)
        Comment = 'Remove DC offset: [All file]';
    elseif any(abs(time) > 2)
        Comment = sprintf('Remove DC offset: [%1.3fs,%1.3fs]', time(1), time(2));
    else
        Comment = sprintf('Remove DC offset: [%dms,%dms]', round(time(1)*1000), round(time(2)*1000));
    end
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
    % Get options
    if isfield(sProcess.options, 'baseline') && isfield(sProcess.options.baseline, 'Value') && iscell(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value{1})
        BaselineBounds = sProcess.options.baseline.Value{1};
    else
        BaselineBounds = [];
    end
    % Raw files: allow to select a segment of recordings that is currently loaded in memory
    if strcmpi(sInput.FileType, 'raw') && (isempty(BaselineBounds) || ((BaselineBounds(1) < sInput.TimeVector(1)) || (BaselineBounds(2) > sInput.TimeVector(end))))
        % Check for previously saved baseline average
        if isfield(sInput, 'RawMeanBaseline') && ~isempty(sInput.RawMeanBaseline)
            meanBaseline = sInput.RawMeanBaseline;
            BaselineRead = sInput.BaselineRead;
        else
            % Progress bar comment
            bst_progress('text', ['Running process: ' sProcess.Comment ' [reading baseline]']);
            % Get the processed channels
            SensorTypes = sProcess.options.sensortypes.Value;
            % Calculate the average over all the time points for the required
            [meanBaseline, BaselineRead] = GetRawTimeAverage(sInput.FileName, sInput.ChannelFile, SensorTypes, BaselineBounds, length(sInput.TimeVector));
            % Save for the next iteration
            sInput.RawMeanBaseline = meanBaseline;
            sInput.BaselineRead = BaselineRead;
        end
    % Regular case: all the baseline is available in the currently loaded segment
    else
        % Get baseline indices
        if ~isempty(BaselineBounds)
            iBaseline = panel_time('GetTimeIndices', sInput.TimeVector, BaselineBounds);
            if isempty(iBaseline)
                bst_report('Error', sProcess, [], 'Invalid baseline definition.');
                sInput = [];
                return;
            end
        % Use all file
        else
            iBaseline = 1:size(sInput.A,2);
        end
        % Compute baseline mean
        meanBaseline = mean(sInput.A(:,iBaseline,:), 2);
        % Baseline that was used for real
        BaselineRead = sInput.TimeVector(iBaseline([1 end]));
    end
    % Remove baseline
    sInput.A = bst_bsxfun(@minus, sInput.A, meanBaseline);
    % History
    sInput.HistoryComment = sprintf('Remove baseline offset: [%1.2f, %1.2f] ms', BaselineRead * 1000);
end


%% ===== GET RAW TIME AVERAGE ====
function [Avg, TimeBoundsRead] = GetRawTimeAverage(DataFile, ChannelFile, SensorTypes, TimeBounds, BlockSize)
    % Read the sFile structure
    DataMat = in_bst_data(DataFile, 'F', 'Time');
    sFile = DataMat.F;
    % Read the channel file
    ChannelMat = in_bst_channel(ChannelFile);
    % Find time samples we are supposed to read
    if ~isempty(TimeBounds)
        iTime = panel_time('GetTimeIndices', DataMat.Time, TimeBounds);
    else
        iTime = 1:length(DataMat.Time);
    end
    % Only support reading from non-epoched files
    iEpoch = 1;
    % Number of blocks to read
    nBlocks = ceil(length(iTime) ./ BlockSize);
    % Get the selected channels
    iChannels = channel_find(ChannelMat.Channel, SensorTypes);
    % Initialize average matrix
    Avg = zeros(length(iChannels), 1);
    % NOTE: FORCE READING CLEAN DATA (CTF compensators + Previous SSP)
    ImportOptions = db_template('ImportOptions');
    ImportOptions.ImportMode = 'Time';
    ImportOptions.DisplayMessages = 0;
    ImportOptions.UseCtfComp = 1;
    ImportOptions.UseSsp = 1;
    ImportOptions.RemoveBaseline = 'no';
    % Loop and read all the indices
    for i = 1:nBlocks
        % Get time indices for this block
        iTimeBlock = min([(i-1)*BlockSize+1; i*BlockSize], length(iTime));
        % Convert to samples in the file
        SamplesBounds = sFile.prop.samples(1) + iTime(iTimeBlock) - 1;
        % Read block of data
        F = in_fread(sFile, ChannelMat, iEpoch, SamplesBounds, iChannels, ImportOptions);
        % Add time points to current average
        Avg = Avg + sum(F,2) ./ length(iTime);
    end
    % Return time bounds that where read from the file
    TimeBoundsRead = DataMat.Time([iTime(1), iTime(end)]);
end




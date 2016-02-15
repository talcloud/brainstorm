function varargout = process_extract_values( varargin )
% PROCESS_EXTRACT_VALUES Extract values from multiple files and concatenate in a matrix file.
% 
% USAGE:   OutputFiles = process_extract_values('Run',     sProcess, sInputs)
%            ConcatMat = process_extract_values('Extract', sProcess, sInputs, OPTIONS)

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
% Authors: Francois Tadel, 2015-2016

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract values';
    sProcess.FileTag     = '| extract';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Extract';
    sProcess.Index       = 350;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/Statistics?highlight=%28Extract+values%29#Histograms';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % === GENERIC EXTRACT OPTIONS
    sProcess = DefineExtractOptions(sProcess);
    
    % === CONCATENATE
    sProcess.options.dim.Comment = {'Concatenate signals (dimension 1)', 'Concatenate time (dimension 2)'};
    sProcess.options.dim.Type    = 'radio';
    sProcess.options.dim.Value   = 2;
    % === OUTPUT COMMENT
    sProcess.options.Comment.Comment = 'Comment (empty=default): ';
    sProcess.options.Comment.Type    = 'text';
    sProcess.options.Comment.Value   = '';
end


%% ===== DEFINE EXTRACT OPTIONS =====
function sProcess = DefineExtractOptions(sProcess)
    % === SELECT: TIME WINDOW
    sProcess.options.timewindow.Comment    = 'Time window:';
    sProcess.options.timewindow.Type       = 'timewindow';
    sProcess.options.timewindow.Value      = [];
    sProcess.options.timewindow.InputTypes = {'data', 'results', 'timefreq', 'matrix'};
    % === SELECT: CHANNELS
    sProcess.options.sensortypes.Comment    = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type       = 'text';
    sProcess.options.sensortypes.Value      = '';
    sProcess.options.sensortypes.InputTypes = {'data'};
    % === SELECT: FREQUENCY RANGE
    sProcess.options.freqrange.Comment    = 'Frequency range: ';
    sProcess.options.freqrange.Type       = 'freqrange';
    sProcess.options.freqrange.Value      = [];
    sProcess.options.freqrange.InputTypes = {'timefreq'};
    % === SELECT: ROWS
    sProcess.options.rows.Comment    = 'Signals names or indices (empty=all): ';
    sProcess.options.rows.Type       = 'text';
    sProcess.options.rows.Value      = '';
    sProcess.options.rows.InputTypes = {'timefreq', 'matrix'};
    
    % === SCOUTS SELECTION ===
    sProcess.options.scoutsel.Comment    = 'Use scouts';
    sProcess.options.scoutsel.Type       = 'scout_confirm';
    sProcess.options.scoutsel.Value      = {};
    sProcess.options.scoutsel.InputTypes = {'results'};
    % === SCOUTS: FUNCTION
    sProcess.options.scoutfunc.Comment    = {'Mean', 'Max', 'PCA', 'Std', 'All', 'Scout function:'};
    sProcess.options.scoutfunc.Type       = 'radio_line';
    sProcess.options.scoutfunc.Value      = 1;
    sProcess.options.scoutfunc.InputTypes = {'results'};
    
    % === NORM XYZ
    sProcess.options.isnorm.Comment    = 'Compute absolute values (or norm for unconstrained sources)';
    sProcess.options.isnorm.Type       = 'checkbox';
    sProcess.options.isnorm.Value      = 0;
    sProcess.options.isnorm.InputTypes = {'results'};
    % === ABSOLUTE VALUE
    sProcess.options.isabs.Comment    = 'Compute absolute values';
    sProcess.options.isabs.Type       = 'checkbox';
    sProcess.options.isabs.Value      = 0;
    sProcess.options.isabs.InputTypes = {'data', 'timefreq', 'matrix'};
    
    % === AVERAGE: TIME
    sProcess.options.avgtime.Comment    = 'Average selected time window';
    sProcess.options.avgtime.Type       = 'checkbox';
    sProcess.options.avgtime.Value      = 0;
    sProcess.options.avgtime.InputTypes = {'data', 'results', 'timefreq', 'matrix'};
    % === AVERAGE: CHANNELS
    sProcess.options.avgrow.Comment    = 'Average selected signals';
    sProcess.options.avgrow.Type       = 'checkbox';
    sProcess.options.avgrow.Value      = 0;
    sProcess.options.avgrow.InputTypes = {'data', 'timefreq', 'matrix'};
    % === AVERAGE: FREQUENCY
    sProcess.options.avgfreq.Comment    = 'Average selected frequency band';
    sProcess.options.avgfreq.Type       = 'checkbox';
    sProcess.options.avgfreq.Value      = 0;
    sProcess.options.avgfreq.InputTypes = {'timefreq'};
end


%% ===== GET EXTRACT OPTIONS =====
function OPTIONS = GetExtractOptions(sProcess, sInputs)
    % Time window
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow) && ~isempty(sProcess.options.timewindow.Value) && iscell(sProcess.options.timewindow.Value)
        OPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    else
        OPTIONS.TimeWindow = [];
    end
    % Sensor type
    if ismember(sInputs(1).FileType, {'data'}) && isfield(sProcess.options, 'sensortypes') && ~isempty(sProcess.options.sensortypes) && ~isempty(sProcess.options.sensortypes.Value)
        OPTIONS.SensorTypes = sProcess.options.sensortypes.Value;
    else
        OPTIONS.SensorTypes = [];
    end
    % Row indices
    if ismember(sInputs(1).FileType, {'results', 'timefreq', 'matrix'}) && isfield(sProcess.options, 'rows') && ~isempty(sProcess.options.rows) && ~isempty(sProcess.options.rows.Value)
        OPTIONS.Rows = sProcess.options.rows.Value;
    else
        OPTIONS.Rows = [];
    end
    % Freq indices
    if ismember(sInputs(1).FileType, {'timefreq'}) && isfield(sProcess.options, 'freqrange') && isfield(sProcess.options.freqrange, 'Value') && iscell(sProcess.options.freqrange.Value) && (length(sProcess.options.freqrange.Value) == 3) && (length(sProcess.options.freqrange.Value{1}) == 2)
        OPTIONS.FreqRange = sProcess.options.freqrange.Value{1};
    else
        OPTIONS.FreqRange = [];
    end
    % Scouts: Selection
    if isfield(sProcess.options, 'scoutsel') && isfield(sProcess.options.scoutsel, 'Value') && isfield(sProcess.options, 'scoutfunc') && isfield(sProcess.options.scoutfunc, 'Value')
        OPTIONS.ScoutSel = sProcess.options.scoutsel.Value;
        switch (sProcess.options.scoutfunc.Value)
            case 1, OPTIONS.ScoutFunc = 'mean';
            case 2, OPTIONS.ScoutFunc = 'pca';
            case 3, OPTIONS.ScoutFunc = 'all';
        end
    else
        OPTIONS.ScoutSel = [];
    end    
    % Absolute values / Norm
    OPTIONS.isAbsolute = 0;
    if isfield(sProcess.options, 'isabs') && isfield(sProcess.options.isabs, 'Value')
        OPTIONS.isAbsolute = sProcess.options.isabs.Value;
    end
    if isfield(sProcess.options, 'isnorm') && isfield(sProcess.options.isnorm, 'Value')
        OPTIONS.isAbsolute = sProcess.options.isnorm.Value;
    end
    % Averages
    if isfield(sProcess.options, 'avgtime') && isfield(sProcess.options.avgtime, 'Value')
        OPTIONS.isAvgTime = sProcess.options.avgtime.Value;
    else
        OPTIONS.isAvgTime = 0;
    end
    if isfield(sProcess.options, 'avgrow') && isfield(sProcess.options.avgrow, 'Value')
        OPTIONS.isAvgRow = sProcess.options.avgrow.Value;
    else
        OPTIONS.isAvgRow = 0;
    end
    if isfield(sProcess.options, 'avgfreq') && isfield(sProcess.options.avgfreq, 'Value')
        OPTIONS.isAvgFreq = sProcess.options.avgfreq.Value;
    else
        OPTIONS.isAvgFreq = 0;
    end
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = [sProcess.Comment, ': [', process_extract_time('GetTimeString', sProcess) ']'];
    if isfield(sProcess.options, 'freqrange') && isfield(sProcess.options.freqrange, 'Value') && iscell(sProcess.options.freqrange.Value) && (length(sProcess.options.freqrange.Value) == 3) && (length(sProcess.options.freqrange.Value{1}) == 2)
        FreqRange = sProcess.options.freqrange.Value{1};
        if (FreqRange(1) == FreqRange(2))
            Comment = [Comment, ' ' num2str(FreqRange(1)) 'Hz'];
        else
            Comment = [Comment, ' ' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) 'Hz'];
        end
    end
    if isfield(sProcess.options, 'sensortypes') && isfield(sProcess.options.sensortypes, 'Value') && ~isempty(sProcess.options.sensortypes.Value)
        Comment = [Comment, ' ', sProcess.options.sensortypes.Value];
    end
    if isfield(sProcess.options, 'rows') && isfield(sProcess.options.rows, 'Value') && ~isempty(sProcess.options.rows.Value)
        Comment = [Comment, ' ', sProcess.options.rows.Value];
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned variable
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    % Get generic options
    OPTIONS = GetExtractOptions(sProcess, sInputs);
    % Output options
    OPTIONS.Dim = sProcess.options.dim.Value;
    
    % ===== EXTRACT =====
    [newMat, FileType] = Extract(sProcess, sInputs, OPTIONS);
    if isempty(newMat)
        return;
    end

    % ===== SAVE FILE =====
    % Get output study
    [sStudy, iStudy, Comment] = bst_process('GetOutputStudy', sProcess, sInputs);
    % Comment: forced in the options
    if isfield(sProcess.options, 'Comment') && isfield(sProcess.options.Comment, 'Value') && ~isempty(sProcess.options.Comment.Value)
        newMat.Comment = sProcess.options.Comment.Value;
    % Comment: Process default
    else
        newMat.Comment = FormatComment(sProcess);
        if (length(sInputs) > 1)
            newMat.Comment = [newMat.Comment, sprintf(' (%d files)', length(sInputs))];
        end
    end
    % Timefreq: the file tag might be longer than the file type
    if strcmpi(FileType, 'timefreq')
        fileTag = bst_process('GetFileTag', sInputs(1).FileName);
    else
        fileTag = FileType;
    end
    % Output filename
    OutputFiles{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), [fileTag '_concat']);
    % Save on disk
    bst_save(OutputFiles{1}, newMat, 'v6');
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, newMat);
end


%% ===== EXTRACT =====
function [newMat, newFileType] = Extract(sProcess, sInputs, OPTIONS)

    % ===== GET INPUTS =====
    % Intialize returned variables
    newMat = [];
    newFileType = [];
    OutValue = cell(1, length(sInputs));
    OutNames = {};
    OutEvents = [];
    ChannelFile = [];
    % Default options
    defOPTIONS = struct(...
        'TimeWindow',  [], ...
        'SensorTypes', [], ...
        'Rows',        [], ...
        'FreqRange',   [], ...
        'ScoutSel',    [], ...
        'ScoutFunc',   [], ...
        'isAbsolute',  0, ...
        'isAvgTime',   0, ...
        'isAvgRow',    0, ...
        'isAvgFreq',   0, ...
        'Dim',         2);      % Dimension on which to concatenate (0,1,2,3,4)
    OPTIONS = struct_copy_fields(OPTIONS, defOPTIONS, 0);
    
    % ===== CHECK INPUTS =====
    % Types of files in input
    inFileType = sInputs(1).FileType;
    % Is time modified with these options
    isTimeModified = ~OPTIONS.isAvgTime && ((OPTIONS.Dim == 1) || (length(sInputs) == 1));
    % Forbid concatenating in dimension 1 multiple source files
    if strcmpi(inFileType, 'results') && (OPTIONS.Dim == 1) && (length(sInputs) > 1) && isempty(OPTIONS.ScoutSel)
        bst_report('Error', sProcess, sInputs(1), 'Cannot concatenate full source files in dimension 1.');
        return;
    end
    % Do not accept connectivity or PAC files
    if ~isempty(strfind(sInputs(1).FileName, '_connect1')) || ~isempty(strfind(sInputs(1).FileName, '_connectn')) || ~isempty(strfind(sInputs(1).FileName, '_pac'))|| ~isempty(strfind(sInputs(1).FileName, '_dpac'))
        bst_report('Error', sProcess, sInputs(1), 'Connectivity and PAC files are not supported yet. Ask on the forum if you need it.');
        return;
    end
    
    % ===== PREPARE HISTORY =====
    % History
    strHistory = ['Extract data: ' FormatComment(sProcess)];
    if OPTIONS.isAvgRow
        strHistory = [strHistory ' AvgRow'];
    end
    if OPTIONS.isAvgTime
        strHistory = [strHistory ' AvgTime'];
    end
    if OPTIONS.isAvgFreq
        strHistory = [strHistory ' AvgFreq'];
    end
    if (length(sInputs) > 1)
        strHistory = [strHistory sprintf(' (%d files)', length(sInputs))];
    end
    % History: Add all the file names
    historyMat.History = [];
    historyMat = bst_history('add', historyMat, 'extract', strHistory);
    if ismember(OPTIONS.Dim, [1 2 3 4])
        historyMat.History = repmat(historyMat.History, length(sInputs) + 1, 1);
        historyMat.History(2:end,3) = cellfun(@(c)cat(2, ' - File: ', c), {sInputs.FileName}, 'UniformOutput', 0);
    end
    
    % ===== LOOP ON FILES =====
    for iInput = 1:length(sInputs)
        bst_progress('text', sprintf('Reading input files... [%d/%d]', iInput, length(sInputs)));
        
        % === LOAD FILE ===
        isAbsApplied = 0;
        % Extract scout
        if strcmpi(sInputs(1).FileType, 'results') && ~isempty(OPTIONS.ScoutSel)
            % Options for LoadInputFile()
            LoadOptions.LoadFull    = 1;    % Load full source results
            LoadOptions.IgnoreBad   = 1;    % Do not read bad segments
            LoadOptions.ProcessName = func2str(sProcess.Function);
            LoadOptions.TargetFunc  = OPTIONS.ScoutFunc;
            LoadOptions.isNorm      = OPTIONS.isAbsolute;
            % Load reference signal
            sLoaded = bst_process('LoadInputFile', sInputs(iInput).FileName, OPTIONS.ScoutSel, OPTIONS.TimeWindow, LoadOptions);
            if isempty(sLoaded) || isempty(sLoaded.Data)
                bst_report('Error', sProcess, sInputs(iInput), 'Could not extract scouts.');
                return;
            end
            % Add components labels
            if (sLoaded.nComponents == 3)
                sLoaded.RowNames = sLoaded.RowNames(:)';
                sLoaded.RowNames = [cellfun(@(c)cat(2,c,'.1'), sLoaded.RowNames, 'UniformOutput', 0); ...
                                    cellfun(@(c)cat(2,c,'.2'), sLoaded.RowNames, 'UniformOutput', 0); ...
                                    cellfun(@(c)cat(2,c,'.3'), sLoaded.RowNames, 'UniformOutput', 0)];
                sLoaded.RowNames = sLoaded.RowNames(:);
            end
            % Convert to matrix structure
            FileMat = db_template('matrixmat');
            FileMat.Value       = sLoaded.Data;
            FileMat.Description = sLoaded.RowNames;
            FileMat.Comment     = sLoaded.Comment;
            FileMat.Time        = sLoaded.Time;
            FileMat.nAvg        = sLoaded.nAvg;
            FileMat.SurfaceFile = sLoaded.SurfaceFile;
            FileMat.Atlas       = sLoaded.Atlas;
            % Interpret as amtrix file in the rest of the function
            inFileType = 'matrix';
            MatValues  = FileMat.Value;
            %isAbsApplied = 1;
        % Regular file
        else
            % Load file 
            [FileMat, matName] = in_bst(sInputs(iInput).FileName, OPTIONS.TimeWindow);
            % File could not be loaded
            if isempty(FileMat) || isempty(FileMat.(matName))
                bst_report('Error', sProcess, sInputs(iInput), 'Could not load anything from the input file. Check the requested time window.');
                return;
            end
            % Unconstrained source files: compute norm
            if strcmpi(inFileType, 'results') && OPTIONS.isAbsolute && (FileMat.nComponents > 1)
                FileMat = process_source_flat('Compute', FileMat, 'rms');
                isAbsApplied = 1;
            end
            % Get values to extract
            MatValues = FileMat.(matName);
        end

        % Get the signals descriptions
        switch (inFileType)
            case 'data'
                % Load channel file (only if new one)
                if ~isequal(ChannelFile, sInputs(iInput).ChannelFile)
                    ChannelFile = sInputs(iInput).ChannelFile;
                    ChannelMat = in_bst_channel(ChannelFile);
                end
                % Select sensors
                if ~isempty(OPTIONS.SensorTypes)
                    iChannels = channel_find(ChannelMat.Channel, OPTIONS.SensorTypes);
                    MatValues = MatValues(iChannels,:,:);
                    Description = {ChannelMat.Channel(iChannels).Name}';
                else
                    Description = {ChannelMat.Channel.Name}';
                end
                Freqs = [];
            case 'results'
                Description = [];
                Freqs = [];
            case 'timefreq'
                Description = FileMat.RowNames;
                % Get file frequency vector
                if iscell(FileMat.Freqs)
                    BandBounds = process_tf_bands('GetBounds', FileMat.Freqs);
                    FreqVector = mean(BandBounds,2);
                else
                    FreqVector = FileMat.Freqs;
                end
                % Rounds the frequency vector, to have the same level of precision as the process (3 significant digits)
                FreqVector = round(FreqVector * 1000) / 1000;
                % Keep only selected frequencies
                if ~isempty(OPTIONS.FreqRange) && ~isempty(FreqVector)
                    iFreqs = find((FreqVector >= OPTIONS.FreqRange(1)) & (FreqVector <= OPTIONS.FreqRange(2)));
                    if isempty(iFreqs)
                        bst_report('Error', sProcess, sInputs(iInput), 'Invalid frequency range.');
                        return;
                    end
                    MatValues = MatValues(:,:,iFreqs);
                    % Keep only the selected frequencies
                    if iscell(FileMat.Freqs)
                        Freqs = FileMat.Freqs(iFreqs,:);
                    else
                        Freqs = FileMat.Freqs(iFreqs);
                    end
                    FreqVector = FreqVector(iFreqs);
                else
                    Freqs = FileMat.Freqs;
                end
            case 'matrix'
                Description = FileMat.Description;
                Freqs = [];
        end
        
        % === SELECT ROWS ===
        if ~isempty(OPTIONS.Rows)
            % Try to find row names in the file
            SelRowNames = strtrim(str_split(OPTIONS.Rows, ',;'));
            iRows = find(ismember(Description, SelRowNames));
            % Else, try selecting by index number
            if isempty(iRows)
                iRows = str2num(OPTIONS.Rows);
                iRows((iRows < 1) | (iRows > length(Description))) = [];
            end
            % Else, try selecting by channel type (for time-frequency computed on channels data)
            if isempty(iRows) && strcmpi(inFileType, 'timefreq') && isequal(FileMat.DataType, 'data') && ~isempty(sInputs(iInput).ChannelFile)
                % Load channel file (only if new one)
                if ~isequal(ChannelFile, sInputs(iInput).ChannelFile)
                    ChannelFile = sInputs(iInput).ChannelFile;
                    ChannelMat = in_bst_channel(ChannelFile);
                end
                % Select sensors
                iChanTF = channel_find(ChannelMat.Channel, OPTIONS.Rows);
                % If something was found: try finding the corresponding channel names in the TF file
                if ~isempty(iChanTF)
                    iRows = find(ismember(Description, {ChannelMat.Channel(iChanTF).Name}));
                end
            end
            % Invalid indices selected
            if (any(iRows > size(MatValues,1)))
                bst_report('Error', sProcess, sInputs(iInput), 'Invalid signals selection.');
                return;
            end
            % If there is something to select
            if ~isempty(iRows)
                MatValues = MatValues(iRows,:,:);
                if ~isempty(Description)
                    Description = Description(iRows);
                end
            end
        else
            iRows = [];
        end
        % Error: Nothing selected
        if isempty(MatValues)
            bst_report('Error', sProcess, sInputs(iInput), 'Could not load anything from the input file. Please check the signal selection.');
            return;
        end
        
        % === ABSOLUTE VALUES ===
        if OPTIONS.isAbsolute && ~isAbsApplied
            MatValues = abs(MatValues);
        end
        % === AVERAGE ROWS ===
        if OPTIONS.isAvgRow && (size(MatValues,1) > 1)
            MatValues = mean(MatValues, 1);
            Description = {'AVG'};
        end
        % === AVERAGE TIME ===
        if OPTIONS.isAvgTime && (size(MatValues,2) > 1)
            % Compute average in time
            MatValues = mean(MatValues, 2);
            % If we are concatenating multiple files in time dimension: keep only one time point for each file
            if (OPTIONS.Dim == 2) && (length(sInputs) > 1)
                FileMat.Time = FileMat.Time(1);
            % If we are not concatenating in time: duplicate the matrix in time so we keep a time segment
            else
                MatValues = [MatValues, MatValues];
                FileMat.Time = [FileMat.Time(1), FileMat.Time(end)];
            end
        end
        % === AVERAGE FREQUENCY ===
        if OPTIONS.isAvgFreq && (size(MatValues,3) > 1) && ~isempty(Freqs)
            MatValues = mean(MatValues, 3);
            Freqs = {'AVG', [num2str(FreqVector(1)), ', ' num2str(FreqVector(end))], 'mean'};
        end

        % === STORE SIGNALS ===
        % Add to load of loaded blocks
        OutValue{iInput} = MatValues;
        % Concatenating in dimension 1: Add the file name
        if (OPTIONS.Dim == 1) && ~isempty(Description)
            Description = cellfun(@(c)cat(2, c, ' @ ', sInputs(iInput).FileName), Description, 'UniformOutput', 0);
        % No concatenation: Update the modified signals list
        elseif (OPTIONS.Dim == 0) && ~isempty(Description)
            switch (inFileType)
                case 'timefreq',  FileMat.RowNames    = Description;
                case 'matrix',    FileMat.Description = Description;
            end
        end
        % First signal
        if (iInput == 1)
            % Save signal names
            OutNames = Description(:);
            % Get time definition
            if (length(FileMat.Time) == 1)
                if (OPTIONS.Dim == 2) && (length(sInputs) > 1)
                    sfreq = 1;
                    tstart = 1;
                else
                    sfreq = 1000;
                    tstart = FileMat.Time(1);
                end
            else 
                sfreq = 1/(FileMat.Time(2) - FileMat.Time(1));
                tstart = FileMat.Time(1);
            end
            % Save first file in memory: will be used as a base (or store all the files if not need of concatenation)
            if (OPTIONS.Dim == 0)
                LoadedMat = cell(1, length(sInputs));
            else
                LoadedMat = cell(1,1);
            end
            LoadedMat{1} = FileMat;
        else
            % Make sure the dimensions are ok
            switch (OPTIONS.Dim)
                case 1
                    isDimOk = (size(OutValue{1},2) == size(OutValue{iInput},2)) && (size(OutValue{1},3) == size(OutValue{iInput},3));
                    % Concatenate names of signals
                    OutNames = cat(1, OutNames, Description(:));
                case 2
                    isDimOk = (size(OutValue{1},1) == size(OutValue{iInput},1)) && (size(OutValue{1},3) == size(OutValue{iInput},3));
                case 3
                    error('Dimension 3 is not supported for concatenation.');
                case 4
                    isDimOk = (size(OutValue{1},1) == size(OutValue{iInput},1)) && (size(OutValue{1},2) == size(OutValue{iInput},2)) && (size(OutValue{1},3) == size(OutValue{iInput},3));
                case 0   % No check, just stacking all the loaded files
                    isDimOk = 1;
                    LoadedMat{iInput} = FileMat;
                otherwise
                    error('Concatenation dimension is not specified.');
            end
            % Error management
            if ~isDimOk
                bst_report('Error', sProcess, sInputs(iInput), sprintf('The dimensions of file #%d do not match file #1.', iInput));
                return;
            end
        end
        
        % === EVENTS ===
        if isfield(FileMat, 'Events') && ~isempty(FileMat.Events) && isTimeModified
            if isempty(OutEvents)
                OutEvents = FileMat.Events;
            else
                sFile.events = OutEvents;
                sFile.prop.sfreq = sfreq;
                sFile = import_events(sFile, [], FileMat.Events);
                OutEvents = sFile.events;
            end
        end
    end
    
    % ===== CONCATENATE SIGNALS ====
    bst_progress('text', 'Concatenating files...');
    % Concatenate data matrix
    if (OPTIONS.Dim ~= 0)
        OutValue = {cat(OPTIONS.Dim, OutValue{:})};
    end
    
    % ===== CREATE STRUCTURE =====
    % Determine the output file format
    switch (inFileType)
        case 'data'
            if OPTIONS.isAvgRow || ((OPTIONS.Dim == 1) && (length(sInputs) > 1)) || (length(ChannelMat.Channel) ~= size(OutValue{1},1))
                newFileType = 'matrix';
            else
                newFileType = inFileType;
            end
        case 'results'
            if OPTIONS.isAvgRow || ((OPTIONS.Dim == 1) && (length(sInputs) > 1)) || ~isempty(iRows)
                newFileType = 'matrix';
            else
                newFileType = inFileType;
            end
        case 'timefreq'
            newFileType = inFileType;
        case 'matrix'
            newFileType = inFileType;
    end
    % If the input and output file types are the same: keep the loaded structure
    if isequal(newFileType, inFileType)
        newMat = LoadedMat;
    % Else: Create a new empty matrix structure
    else
        newMat = {db_template('matrixmat')};
        if (OPTIONS.Dim == 0)
            newMat = repmat(newMat, 1, length(LoadedMat));
        end
    end
    % Process all the structures in the same way
    for i = 1:length(LoadedMat)
        % Replicate value in time for display if necessary
        if ismember(newFileType, {'data', 'results', 'matrix'}) && (size(OutValue{i},2) == 1)
            OutValue{i} = [OutValue{i}, OutValue{i}];
        end
        % Copy the values in the appropriate fields
        switch (newFileType)
            case 'data'
                newMat{i}.F      = OutValue{i};
                newMat{i}.Events = OutEvents;
            case 'results'
                newMat{i}.ImageGridAmp = OutValue{i};
                % Disconnect from parent file if the time was altered
                if isTimeModified || (length(sInputs) > 1)
                    newMat{i}.DataFile = [];
                end
            case 'timefreq'
                newMat{i}.TF       = OutValue{i};
                newMat{i}.Freqs    = Freqs;
                % If there is not concatenation: do not force the row names to be all the same
                if (OPTIONS.Dim ~= 0)
                    newMat{i}.RowNames = OutNames;
                end
                % Disconnect from parent file if the time was altered
                if (length(sInputs) > 1)
                    newMat{i}.DataFile = [];
                end
                % Change the type of data if necessary
                switch (newMat{i}.DataType)
                    case 'data'
                        if OPTIONS.isAvgRow || ((OPTIONS.Dim == 1) && (length(sInputs) > 1)) || ~isempty(iRows)
                            newMat{i}.DataType = 'matrix';
                        end
                    case 'results'
                        if OPTIONS.isAvgRow || ((OPTIONS.Dim == 1) && (length(sInputs) > 1)) || ~isempty(iRows)
                            newMat{i}.DataType = 'matrix';
                        end
                end
            case 'matrix'
                newMat{i}.Value       = OutValue{i};
                newMat{i}.Events      = OutEvents;
                % If there is not concatenation: do not force the row names to be all the same
                if (OPTIONS.Dim ~= 0)
                    newMat{i}.Description = OutNames;
                end
        end
        % Final time vector
        newMat{i}.Time = (0:size(OutValue{i},2)-1) ./ sfreq + tstart;
        % Add common history field
        newMat{i} = bst_history('add', newMat{i}, historyMat.History);
    end
    % Just return one structure
    if (OPTIONS.Dim ~= 0)
        newMat = newMat{1};
    end
end



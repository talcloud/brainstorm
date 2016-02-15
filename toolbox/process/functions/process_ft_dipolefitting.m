function varargout = process_ft_dipolefitting( varargin )
% PROCESS_FT_DIPOLEFITTING: Call FieldTrip function ft_dipolefitting.
%
% REFERENCES: 
%     - http://www.fieldtriptoolbox.org/reference/ft_dipolefitting
%     - http://www.fieldtriptoolbox.org/tutorial/natmeg/dipolefitting

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
% Authors: Jeremy T. Moreau, Elizabeth Bock, Francois Tadel, 2015-2016

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'FieldTrip: ft_dipolefitting';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 350;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/DipoleFitting';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Comment
    sProcess.options.filetag.Comment = 'File tag: ';
    sProcess.options.filetag.Type    = 'text';
    sProcess.options.filetag.Value   = '';
    % Options: Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Options: Sensor type
    sProcess.options.sensortypes.Comment = 'Sensor type or names: ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG';
    % Options: Dipole model
    sProcess.options.dipolemodel.Comment = {'Moving dipole', 'Regional dipole', 'Dipole type: '};
    sProcess.options.dipolemodel.Type    = 'radio_line';
    sProcess.options.dipolemodel.Value   = 1;
    % Options: Number of dipoles
    sProcess.options.numdipoles.Comment = 'Number of dipoles to fit: ';
    sProcess.options.numdipoles.Type    = 'value';
    sProcess.options.numdipoles.Value   = {1, 'dipoles', 0};
    % Options: Symmetry constraint
    sProcess.options.symmetry.Comment = 'Left-right symmetry constraint (for two dipoles only)';
    sProcess.options.symmetry.Type    = 'checkbox';
    sProcess.options.symmetry.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFile = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFile = [];
    % Initialize fieldtrip
    bst_ft_init();
    
    % ===== GET OPTIONS =====
    SensorTypes = sProcess.options.sensortypes.Value;
    TimeWindow  = sProcess.options.timewindow.Value{1};
    NumDipoles  = sProcess.options.numdipoles.Value{1};
    % Dipole model
    switch (sProcess.options.dipolemodel.Value)
        case 1,  DipoleModel = 'moving';
        case 2,  DipoleModel = 'regional';
    end
    % Symmetry constraints
    if (sProcess.options.symmetry.Value) && (NumDipoles == 2)
        SymmetryConstraint = 'y';
    else
        SymmetryConstraint = [];
    end
    % File tag
    FileTag = sProcess.options.filetag.Value;
    if isempty(FileTag)
        c = clock();
        FileTag = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
    end

    % ===== LOAD INPUTS =====
    % Load channel file
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    % Get selected sensors
    iChannels = channel_find(ChannelMat.Channel, SensorTypes);
    if isempty(iChannels)
        bst_report('Error', sProcess, sInput, ['Channels "' SensorTypes '" not found in channel file.']);
        return;
    end
    % Check the sensor types (only one type allowed)
    AllTypes = unique({ChannelMat.Channel(iChannels).Type});
    if (length(AllTypes) > 1) && all(ismember(AllTypes, {'MEG MAG', 'MEG GRAD'}))
        AllTypes = setdiff(AllTypes, {'MEG MAG', 'MEG GRAD'});
        AllTypes = union(AllTypes, 'MEG');
    end
    if (length(AllTypes) ~= 1)
        bst_report('Error', sProcess, sInput, 'FieldTrip dipole fitting works only on one sensor type at a time.');
        return;
    elseif ~ismember(AllTypes{1}, {'MEG','EEG','MEG MAG','MEG GRAD', 'MEG GRAD2', 'MEG GRAD3'})
        bst_report('Error', sProcess, sInput, 'Only MEG and EEG sensor types are supported at this moment.');
        return;
    end

    % Get head model
    sHeadModel = bst_get('HeadModelForStudy', sInput.iStudy);
    % Error: No head model
    if isempty(sHeadModel)
        bst_report('Error', sProcess, sInput, 'No head model available for this data file.');
        return;
    % Error: Ensure that the correct head model has been selected
    elseif strcmp(SensorTypes,'MEG') && ~(strcmp(sHeadModel.MEGMethod, 'meg_sphere') || strcmp(sHeadModel.MEGMethod, 'os_meg'))
        bst_report('Error', sProcess, sInput, 'Only single sphere or overlapping spheres head models are supported for MEG dipole fitting.');
        return;
    elseif strcmp(SensorTypes,'EEG') && ~strcmp(sHeadModel.EEGMethod, 'eeg_3sphereberg')
        bst_report('Error', sProcess, sInput, 'Only 3-shell spheres head models are supported for EEG dipole fitting.');
        return;
    end
    
    % Load data
    DataFile = sInput.FileName;
    DataMat = in_bst_data(DataFile);
    % Remove bad channels
    iBadChan = find(DataMat.ChannelFlag == -1);
    iChannels = setdiff(iChannels, iBadChan);
    % Error: All channels tagged as bad
    if isempty(iChannels)
        bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
        return;
    end
    

    % ===== CALL FIELDTRIP =====
    % Convert head model and data
    [ftHeadModel, HeadModelMat] = out_fieldtrip_headmodel(sHeadModel.FileName, ChannelMat, iChannels);
    ftData = out_fieldtrip_data(DataFile, ChannelMat, iChannels, 1);
    % Generate rough grid for first estimation
    GridOptions.Method        = 'isotropic';
    GridOptions.nVerticesInit = 4000;
    GridOptions.Resolution    = 0.010;
    GridLoc = bst_sourcegrid(GridOptions, HeadModelMat.SurfaceFile);
    
    % Initialise unlimited progress bar
    bst_progress('start', 'ft_dipolefitting', 'Calling FieldTrip function: ft_dipolefitting...');
    % Prepare FieldTrip cfg structure
    cfg = [];
    cfg.channel     = {ChannelMat.Channel(iChannels).Name};
    cfg.headmodel   = ftHeadModel;
    cfg.latency     = TimeWindow;
    cfg.numdipoles  = NumDipoles;
    cfg.model       = DipoleModel;
    cfg.nonlinear   = 'yes';
    cfg.gridsearch  = 'yes';
    cfg.grid.pos    = GridLoc;
    cfg.grid.inside = ones(size(GridLoc,1),1);
    cfg.grid.unit   = 'm';
    cfg.symmetry    = SymmetryConstraint;
    cfg.feedback    = 'textbar';
    % Optimization function
    if exist('fminunc', 'file')
        cfg.dipfit.optimfun = 'fminunc';
    else 
        cfg.dipfit.optimfun = 'fminsearch';
    end
    % Run ft_dipolefitting
    ftDipole = ft_dipolefitting(cfg, ftData);
    % Check if something was returned
    if isempty(ftDipole) || isempty(ftDipole.dip) || all(ftDipole.dip(1).pos(:) == 0) || all(ftDipole.dip(1).mom(:) == 0)
        bst_report('Error', sProcess, sInput, 'Something went wrong during the execution of the FieldTrip function. Check the command window...');
        return;
    end
    % Get output dipoles
    nTime = length(ftDipole.time);
    dipTime = ftDipole.time;
    switch (DipoleModel)
        case 'moving'
            dipPos = cat(1, ftDipole.dip.pos)';
            dipMom = reshape(cat(2, ftDipole.dip.mom), 3, []);
            dipRv  = cat(2, ftDipole.dip.rv);
        case 'regional'
            dipPos = repmat(ftDipole.dip.pos, nTime, 1)';
            dipMom = reshape(ftDipole.dip.mom, 3, []);
            dipRv  = ftDipole.dip.rv;
    end
    % Replace single values for multiple dipoles
    if (NumDipoles > 1)
        dipRv = reshape(repmat(dipRv, NumDipoles, 1), 1, []);
        dipTime = reshape(repmat(dipTime, NumDipoles, 1), 1, []);        
    end
    
    % ===== OUTPUT STRUCTURE =====
    % Intialize dipole structure
    DipolesMat = db_template('dipolemat');
    DipolesMat.Time        = ftDipole.time;
    DipolesMat.DataFile    = sInput.FileName;
    DipolesMat.DipoleNames = {'ft_dipolefitting'};
    DipolesMat.Subset      = 1;
    DipolesMat.Comment = sprintf('ft_dipolefitting [%ddip,%dtime]: %s', NumDipoles, nTime, FileTag);
    % Estimated dipoles
    for iDip = 1:length(dipTime)
        DipolesMat.Dipole(iDip).Index     = 1;
        DipolesMat.Dipole(iDip).Time      = dipTime(iDip);
        DipolesMat.Dipole(iDip).Origin    = [0, 0, 0];
        DipolesMat.Dipole(iDip).Loc       = dipPos(:,iDip);
        DipolesMat.Dipole(iDip).Amplitude = dipMom(:,iDip);
        DipolesMat.Dipole(iDip).Errors    = 0;
        DipolesMat.Dipole(iDip).Goodness  = 1 - dipRv(iDip);    % RV = RESIDUAL VARIANCE
    end
    % Save FieldTrip configuration structure
    DipolesMat.cfg = cfg;
    if isfield(DipolesMat.cfg, 'headmodel')
        DipolesMat.cfg = rmfield(DipolesMat.cfg, 'headmodel');
    end
    if isfield(ftDipole, 'cfg') && isfield(ftDipole.cfg, 'version')
        DipolesMat.cfg.version = ftDipole.cfg.version;
    end
       
    
    % ===== SAVE FILE =====
    % Create a file with the same name as the input data file
    [fPath, fBase, fExt] = bst_fileparts(file_fullpath(sInput.FileName));
    DipoleFile = file_unique(bst_fullfile(fPath, [strrep(fBase, 'data_', 'dipoles_'), '.mat']));
    % Save new file in Brainstorm format
    bst_save(DipoleFile, DipolesMat);
    
    % ===== UPDATE DATABASE =====
    % Create structure
    BstDipolesMat = db_template('Dipoles');
    BstDipolesMat.FileName = file_short(DipoleFile);
    BstDipolesMat.Comment  = DipolesMat.Comment;
    BstDipolesMat.DataFile = DipolesMat.DataFile;
    % Add to study
    sStudy = bst_get('Study', sInput.iStudy);
    iDipole = length(sStudy.Dipoles) + 1;
    sStudy.Dipoles(iDipole) = BstDipolesMat;
    % Save study
    bst_set('Study', sInput.iStudy, sStudy);
    % Update tree
    panel_protocols('UpdateNode', 'Study', sInput.iStudy);
    % Save database
    db_save();
    
    % Return the input file (as we cannot handle the dipole files in the pipeline editor)
    OutputFile = sInput.FileName;
end




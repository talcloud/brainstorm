function varargout = process_dipole_scanning( varargin )
% PROCESS_DIPOLE_SCANNING: Generates a brainstorm dipole file from the GLS and GLS-P inverse solutions.

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
% Authors: Elizabeth Bock, John C. Mosher, Francois Tadel, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Dipole scanning';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 327;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/TutDipScan';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'dipoles'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
     % === Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
%     % === fit frequency
%     sProcess.options.downsample.Comment = 'Time between dipoles: ';
%     sProcess.options.downsample.Type    = 'value';
%     sProcess.options.downsample.Value   = {0.000, 'ms', []};
    % === Separator
    sProcess.options.sep2.Type = 'separator';
    sProcess.options.sep2.Comment = ' ';
    % === CLUSTERS
    sProcess.options.scouts.Comment = 'Limit scanning to selected scouts';
    sProcess.options.scouts.Type    = 'scout_confirm';
    sProcess.options.scouts.Value   = {};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};
    % === Get options
     % Time window to process
     %FitPeriod  = sProcess.options.downsample.Value{1};
    TimeWindow = sProcess.options.timewindow.Value{1};
    FitPeriod  = 0;
    AtlasList  = sProcess.options.scouts.Value;

    % === Get the sources
    % Read the source file
    sResultP = in_bst_results(sInput.FileName, 0);
    % Get the scouts structures
    if ~isempty(AtlasList)
        [sScouts, AtlasNames, sSurf] = process_extract_scout('GetScoutsInfo', sProcess, sInput, sResultP.SurfaceFile, AtlasList);
    else
        sScouts = [];
    end

    % === Get the results
    if ~isempty(strfind(sResultP.Function, 'gls_p')) || ~isempty(strfind(sResultP.Function, 'glsp'))
       ScanType = 'GLSP';    % works for all three cases of Mosher
    elseif ~isempty(strfind(sResultP.Function, 'mnej')) || ~isempty(strfind(sResultP.Function, 'mnp')) 
       ScanType = 'MNEJP';
    elseif ~isempty(strfind(sResultP.Function, 'glsr'))
       ScanType = 'GLSRP';
    elseif ~isempty(strfind(sResultP.Function, 'dspm'))
       ScanType = 'dSPM';
    elseif ~isempty(strfind(sResultP.Function, 'sloreta'))
       ScanType = 'sLORETA';
    elseif ~isempty(strfind(sInput.FileName, 'zscore')) || ~isempty(strfind(sResultP.Function, 'zscore'))
       ScanType = 'zscore';
    else
       bst_report('Error', sProcess, [], sprintf('%s\n%s', 'Dipole scanning is only available on a performance image matrix', 'i.e.  GLSP, dSPM, sLORETA and MN zscore'));
       return;
    end
    
    if isempty(sResultP.DataFile)
        DataMatP.Time = sResultP.Time;
        DataMatP.F = [];
        DataMatP.Comment = sResultP.Comment;
    else
        DataMatP = in_bst_data(sResultP.DataFile);
    end
    
    if isempty(sResultP.ImageGridAmp)
        sResultP.ImageGridAmp = sResultP.ImagingKernel * DataMatP.F(sResultP.GoodChannel,:); 
    end
    
    % === Get the time
    if ~isempty(TimeWindow)
        SamplesBounds = bst_closest(TimeWindow, DataMatP.Time);
    else
        SamplesBounds = [1, size(sResultP.ImageGridAmp,2)];
    end
    timeVector = DataMatP.Time(SamplesBounds(1):SamplesBounds(2));
    
    % The "performance" image matrix
    P = sResultP.ImageGridAmp(:,SamplesBounds(1):SamplesBounds(2));
    % For contrained or mixed models this will flatten to norm
    sResultPFlat = process_source_flat('Compute',sResultP, 'rms');
    Pscan = sResultPFlat.ImageGridAmp(:,SamplesBounds(1):SamplesBounds(2));
    
    % === Find the index of 'best fit' at every time point
    % Get the selected scouts
    if ~isempty(sScouts)
        scoutVerts = [];
        for iScout = 1:length(sScouts)
            scoutVerts = [scoutVerts sScouts(iScout).Vertices];
        end
        [mag,ind] = max(Pscan(scoutVerts,:));
        maxInd = scoutVerts(ind);
    else
        % don't use scouts, use all vertices
        [mag,maxInd] = max(Pscan,[],1);
    end
    
    % === Prepare a mask for downsampling the number of dipoles to save
    NumDipoles = size(P,2);
    if FitPeriod > 0
        tTotal = timeVector(end) - timeVector(1);
        nNewDipoles = tTotal/FitPeriod;
        dsFactor = round(NumDipoles/nNewDipoles);
    else
        dsFactor = 1;
    end
    temp = zeros(1,dsFactor);
    temp(1) = 1;
    dsMask = repmat(temp, 1, floor(NumDipoles/dsFactor));
    dsMask = logical([dsMask zeros(1,NumDipoles-length(dsMask))]);
    % downsample 
    dsTimeVector = timeVector(dsMask);
    dsMaxInd = maxInd(dsMask);

    % Get headmodel
    sHeadModel = bst_get('HeadModelForStudy', sInput.iStudy);
    if isempty(sHeadModel)
        HeadModelFile = sResultP.HeadModelFile; 
        if isempty(HeadModelFile)
            error('No headmodel available for this result file.');
        end
    else
        HeadModelFile = sHeadModel.FileName;
    end
    HeadModelMat = in_headmodel_bst(HeadModelFile, 0, 'Gain', 'GridLoc', 'GridOrient');
    
    % === find the orientations
    switch (sResultP.nComponents)
        case 0  % MIXED HEAD MODEL
            %fullOrient = zeros(3, size(P,1), size(P,2));
            nTime = size(P,2);
            fullOrient = {};
            % Loop on all the regions
            for iScout = 1:length(sResultP.GridAtlas.Scouts)
                % Get the grid indices for this scouts
                iGrid = sResultP.GridAtlas.Scouts(iScout).GridRows;
                % If no vertices to read from this region: skip
                if isempty(iGrid)
                    continue;
                end
                % Get correpsonding row indices based on the type of region (constrained or unconstrained)
                switch (sResultP.GridAtlas.Scouts(iScout).Region(3))
                    case 'C'
                        fullOrient{end+1} = repmat(sResultP.GridOrient(iGrid,:)', 1, 1, nTime);
                    case {'U','L'}
                        % Convert from indices in GridLoc to indices in the source matrix
                        iSourceRows = bst_convert_indices(iGrid, sResultP.nComponents, sResultP.GridAtlas, 0);
                        fullOrient{end+1} = reshape(P(iSourceRows,:), 3, [], size(P,2));
                    otherwise
                        error(['Invalid region "' sResultP.GridAtlas.Scouts(iScout).Region '"']);
                end
            end
            % Concatenate the blocks
            fullOrient = cat(2, fullOrient{:});
            % Dipole orientation = amplitude in the three orientations
            for jj = 1:length(dsMaxInd)
                orient(:,jj) = fullOrient(:,dsMaxInd(jj),jj);
            end
        case 1
            % use the headmodel orientations
            fullOrient = HeadModelMat.GridOrient';
            for jj = 1:length(dsMaxInd)
                orient(:,jj) = fullOrient(:,dsMaxInd(jj));
            end
        case 2
            error('Not supported.');
        case 3
            % use the dipole orientations
            fullOrient = reshape(P,sResultP.nComponents,[],size(P,2));
            for jj = 1:length(dsMaxInd)
                orient(:,jj) = fullOrient(:,dsMaxInd(jj),jj);
            end
    end

    % === Performace measures for GLSP
    if (strcmp(ScanType, 'GLSP') || strcmp(ScanType,'MNEJP') || strcmp(ScanType,'GLSRP') || strcmp(ScanType,'zscore')) && ~isempty(DataMatP.F)
        isComputePerformanceMeasures = 1;
    end
    
    isComputePerformanceMeasures = 0; % TODO waiting for John to fix performance measures.
    if isComputePerformanceMeasures 
        % === Goodness of fit
        % square the performance at every source, therefore resulting in a
        % scalar squared performance value at every dipolar source
        %P2 = sum(reshape(abs(P).^2,sResultP.nComponents,[]),1);
        %P2 = reshape(P2,[],size(P,2));
        P2 = Pscan .^ 2;

        % get the squared norm of the whitened data
        wd2 = sum(abs(sResultP.Whitener * DataMatP.F(sResultP.GoodChannel,SamplesBounds(1):SamplesBounds(2))).^2,1);

        % the goodness of fit is now calculated by dividing the norm into each
        % performance value
        gof = P2 * diag(1./wd2);

        % === Chi-square
        % the chi square is the difference of the norm and the performance
        % resulting in the error for every source at every time point
        ChiSquare = repmat(wd2,size(P2,1),1) - P2;

        % The reduced chi-square is found by dividing by the degrees of freedom in
        % the error, which (for now) is simply a scalar, since we assume all
        % sources have the same degrees of freedom. Thus ROI modeling will require
        % that all ROIs have the same DOF. 
        DOF = size(sResultP.ImagingKernel,2) - sResultP.nComponents;
        RChiSquare = ChiSquare / DOF;
        
        % downsample
        dsChiSquare = ChiSquare(:,dsMask);
        dsRChiSquare = RChiSquare(:,dsMask);
        dsGOF = gof(:,dsMask);
    
    else
       dsChiSquare = [];
       dsRChiSquare = [];
       dsGOF = [];
       DOF = [];
    end
    
    dsP = Pscan(:,dsMask);

    
    % === Error mask
%     ERR_THRESH = 3; 
%     Error_Mask = RChiSquare < ERR_THRESH;
% 
%     P2mask = zeros(size(P2));
%     P2mask(Error_Mask) = P2(Error_Mask);
% 
%     GOFmask = zeros(size(GOF));
%     GOFmask(Error_Mask) = GOF(Error_Mask);
    
    NumDipoles = length(dsMaxInd);
    
    %% === CREATE OUTPUT STRUCTURE ===
    bst_progress('start', 'Dipole File', 'Saving result...');
    % Get output study
    [sStudy, iStudy] = bst_get('Study', sInput.iStudy);
    % Comment: forced in the options
    if isfield(sProcess.options, 'Comment') && isfield(sProcess.options.Comment, 'Value') && ~isempty(sProcess.options.Comment.Value)
        Comment = sProcess.options.Comment.Value;
    % Comment: process default
    else
        Comment = [DataMatP.Comment ' | ' ScanType '-dipole-scan'];
    end
    % Get base filename
    [fPath, fBase, fExt] = bst_fileparts(sInput.FileName);
    % Create base structure
    DipolesMat = db_template('dipolemat');
    DipolesMat.Comment = Comment;
    DipolesMat.Time    = unique(dsTimeVector);
    
    % Fill structure    
    for i = 1:NumDipoles
        DipolesMat.Dipole(i).Index          = 1;
        DipolesMat.Dipole(i).Time           = dsTimeVector(i);
        DipolesMat.Dipole(i).Origin         = [0 0 0];
        DipolesMat.Dipole(i).Loc            = HeadModelMat.GridLoc(dsMaxInd(i),:)';
        DipolesMat.Dipole(i).Amplitude      = orient(:,i);
        DipolesMat.Dipole(i).Errors         = 0;
        DipolesMat.Dipole(i).Noise          = [];
        DipolesMat.Dipole(i).SingleError    = [];
        DipolesMat.Dipole(i).ErrorMatrix    = [];
        DipolesMat.Dipole(i).ConfVol        = [];
        DipolesMat.Dipole(i).Probability    = [];
        DipolesMat.Dipole(i).NoiseEstimate  = [];
        DipolesMat.Dipole(i).Perform        = dsP(dsMaxInd(i),i);
        
        if ~isempty(dsGOF)
           DipolesMat.Dipole(i).Goodness       = dsGOF(dsMaxInd(i),i);
        end
        if ~isempty(dsChiSquare)
           DipolesMat.Dipole(i).Khi2           = dsChiSquare(dsMaxInd(i),i);
        end
        if ~isempty(DOF)
           DipolesMat.Dipole(i).DOF            = DOF;
        end
    end

    % Create the dipoles names list
    dipolesList = unique([DipolesMat.Dipole.Index]); %unique group names
    DipolesMat.DipoleNames = cell(1,length(dipolesList));
    k = 1; %index of names for groups with subsets
    nChanSet = 1;
    for i = 1:(length(dipolesList)/nChanSet)
        % If more than one channel subset, name the groups according to index
        % and subset number
        if nChanSet > 1
            for j = 1:nChanSet
                DipolesMat.DipoleNames{k} = sprintf('Group #%d (%d)', dipolesList(i), j);
                DipolesMat.Subset(k) = j;
                k=k+1;
            end
        % if only one subsets, name the groups according to index
        else
            DipolesMat.DipoleNames{i} = sprintf('Group #%d', dipolesList(i));
            DipolesMat.Subset(i) = 1;
        end

    end
    DipolesMat.Subset = unique(DipolesMat.Subset);


    % Add data file
    DipolesMat.DataFile = sInput.FileName;
    % Add History field
    DipolesMat = bst_history('add', DipolesMat, 'generate', ['Generated from: ' sInput.FileName]);
   
    %% ===== SAVE NEW FILE =====
    % Get imported base name
    [tmp__, importedBaseName, importedExt] = bst_fileparts(sInput.FileName);
    importedBaseName = strrep(importedBaseName, 'dipoles_', '');
    importedBaseName = strrep(importedBaseName, '_dipoles', '');
    importedBaseName = strrep(importedBaseName, 'dipoles', '');
    % Limit number of chars
    if (length(importedBaseName) > 15)
        importedBaseName = importedBaseName(1:15);
    end
    % Create output filename
    ProtocolInfo = bst_get('ProtocolInfo');
    DipoleFile = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sStudy.FileName), 'dipoles_fit.mat');
    DipoleFile = file_unique(DipoleFile);
    % Save new file in Brainstorm format
    bst_save(DipoleFile, DipolesMat);
    
    
    %% ===== UPDATE DATABASE =====
    % Create structure
    BstDipolesMat = db_template('Dipoles');
    BstDipolesMat.FileName = file_short(DipoleFile);
    BstDipolesMat.Comment  = Comment;
    BstDipolesMat.DataFile = sInput.FileName;
    % Add to study
    iDipole = length(sStudy.Dipoles) + 1;
    sStudy.Dipoles(iDipole) = BstDipolesMat;
    
    % Save study
    bst_set('Study', iStudy, sStudy);
    % Update tree
    panel_protocols('UpdateNode', 'Study', iStudy);
    % Save database
    db_save();

    OutputFiles{1} = DipoleFile;
end






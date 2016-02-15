function varargout = panel_inverse_2015(varargin)
% PANEL_INVERSE_2015: Inverse modeling GUI.
%
% USAGE:  bstPanelNew = panel_inverse_2015('CreatePanel')

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
% Authors: Francois Tadel, 2008-2015

macro_methodcall;
end


%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(Modalities, isShared, HeadModelType) %#ok<DEFNU>
    panelName = 'InverseOptions';
    % Java initializations
    import java.awt.*;
    import javax.swing.*;
    global ALLOW_2015;
    
    % GUI CALL:  panel_wmne('CreatePanel', Modalities, isShared, HeadModelType)
    if iscell(Modalities)
        isProcess = 0;
        isFull = 0;
        % Get default inverse options
        OPTIONS = bst_inverse_linear_2015();
    % PROCESS CALL:  panel_wmne('CreatePanel', sProcess, sFiles)
    else
        isProcess = 1;
        % Get inputs
        sProcess = Modalities;
        sFiles   = isShared;
        % Get inverse options
        OPTIONS = sProcess.options.inverse.Value;
        % List of sensors
        Modalities = intersect(sFiles(1).ChannelTypes, {'MEG MAG', 'MEG GRAD', 'MEG', 'EEG', 'ECOG', 'SEEG'});
        if any(ismember({'MEG MAG','MEG GRAD'}, Modalities))
            Modalities = setdiff(Modalities, 'MEG');
        end
        % Shared kernel
        isShared = (sProcess.options.output.Value == 1);
        isFull   = (sProcess.options.output.Value == 3);
        % Head model type: Get from head model
        sStudyChan = bst_get('ChannelFile', sFiles(1).ChannelFile);
        if ~isempty(sStudyChan) && ~isempty(sStudyChan.HeadModel)
            HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
            HeadModelMat = in_headmodel_bst(HeadModelFile, 0, 'HeadModelType');
            HeadModelType = HeadModelMat.HeadModelType;
        else
            HeadModelType = 'surface';
        end
    end
    % Initializations
    isFirstCombinationWarning = 1;
    
    % ==== FRAME STRUCTURE ====
    % Create main panel: split in top (interface) / left(standard options) / right (details)
    jPanelNew = gui_component('panel');
    jPanelNew.setBorder(BorderFactory.createEmptyBorder(12,12,12,12));
    % Create top panel
    jPanelTop = gui_river([1,1], [0,6,6,6]);
    jPanelNew.add(jPanelTop, BorderLayout.NORTH);
    % Create left panel
    jPanelLeft = java_create('javax.swing.JPanel');
    jPanelLeft.setLayout(BoxLayout(jPanelLeft, BoxLayout.Y_AXIS));
    jPanelNew.add(jPanelLeft, BorderLayout.WEST);
    % Create right panel
    jPanelRight = java_create('javax.swing.JPanel');
    %jPanelRight.setLayout(BoxLayout(jPanelRight, BoxLayout.Y_AXIS));
    jPanelRight.setLayout(GridBagLayout());
    jPanelRight.setBorder(BorderFactory.createEmptyBorder(0,12,0,0));
    jPanelNew.add(jPanelRight, BorderLayout.EAST);
    % Default constrains
    c = GridBagConstraints();
    c.fill    = GridBagConstraints.HORIZONTAL;
    c.weightx = 1;
    c.weighty = 0;
    
    % ==== TOP PANEL ====
    % Warning
    if isempty(ALLOW_2015) || ~ALLOW_2015
        gui_component('label', jPanelTop, [], '<HTML><B><FONT color="#CC0000">Warning: This interface is under development.<BR>Use at your own risk.</B></FONT><BR><BR>', [], [], [], []);
    end
    % Comment
    gui_component('label', jPanelTop, 'br', 'Comment:  ', [], [], [], []);
    jTextComment = gui_component('text', jPanelTop, 'hfill', '', [], '', [], []);
    if isProcess && isfield(OPTIONS, 'Comment') && ~isempty(OPTIONS.Comment)
        jTextComment.setText(OPTIONS.Comment);
    end
    % Linear / Non-linear
    jGroupLinear  = ButtonGroup(); 
    jToogleLinear = gui_component('ToolbarToggle', jPanelTop, 'br', 'Linear', jGroupLinear, '', @(h,ev)UpdatePanel(1), []);
    jToogleNonLin = gui_component('ToolbarToggle', jPanelTop, [], 'Non-linear', jGroupLinear, '', @(h,ev)UpdatePanel(1), []);
    jToogleLinear.setSelected(1);
    % Disable for shared/volume
    if ~strcmpi(HeadModelType, 'surface') || isShared || (exist('isdeployed', 'builtin') && isdeployed)
        jToogleNonLin.setEnabled(0);
    end
    
    % ======================================================================================================
    
    % ==== PANEL: NON-LINEAR ====
    jPanelNonLin = gui_river([1,1], [0,6,6,6], 'Non-linear methods');
        jGroupNonLin  = ButtonGroup();       
        jRadioMem = gui_component('radio', jPanelNonLin, [], 'MEM: Maximum entropy on the mean', jGroupNonLin, '', @(h,ev)UpdatePanel(), []);
        jRadioMem.setSelected(1);
    jPanelLeft.add(jPanelNonLin);
    
    % ==== PANEL: METHOD ====
    jPanelMethod = gui_river([1,1], [0,6,6,6], 'Method');
        jGroupMethod  = ButtonGroup();
        jRadioMethodMn  = gui_component('radio', jPanelMethod, [],   'Minimum norm imaging', jGroupMethod, '', @(h,ev)UpdatePanel(1), []);
        jRadioMethodDip = gui_component('radio', jPanelMethod, 'br', 'Dipole modeling',      jGroupMethod, '', @(h,ev)UpdatePanel(1), []);
        jRadioMethodBf  = gui_component('radio', jPanelMethod, 'br', 'LCMV beamformer',      jGroupMethod, '', @(h,ev)UpdatePanel(1), []);
        % Default selection
        switch lower(OPTIONS.InverseMethod)
            case 'minnorm',  jRadioMethodMn.setSelected(1);
            case 'gls',      jRadioMethodDip.setSelected(1);
            case 'lcmv',     jRadioMethodBf.setSelected(1);
            case 'mem',      disp('BST> Warning: Running MEM from a script is not handled yet.');
        end
    jPanelLeft.add(jPanelMethod);
    
    % ==== PANEL: MEASURE MIN NORM ====
    jPanelMeasureMN = gui_river([1,1], [0,6,6,6], 'Measure');
        jGroupMnMeasure = ButtonGroup();
        jRadioMnCurrent = gui_component('radio', jPanelMeasureMN, [],   'Current density map',  jGroupMnMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioMnDspm    = gui_component('radio', jPanelMeasureMN, 'br', 'dSPM',                 jGroupMnMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioMnSloreta = gui_component('radio', jPanelMeasureMN, 'br', 'sLORETA',              jGroupMnMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioMnPerform = gui_component('radio', jPanelMeasureMN, 'br', 'NP performance index', jGroupMnMeasure, '', @(h,ev)UpdatePanel(), []);
        % Default selection
        switch lower(OPTIONS.InverseMeasure)
            case 'amplitude',    jRadioMnCurrent.setSelected(1);
            case 'performance',  jRadioMnPerform.setSelected(1);
            case 'dspm',         jRadioMnDspm.setSelected(1);
            case 'sloreta',      jRadioMnSloreta.setSelected(1);
            otherwise,           jRadioMnCurrent.setSelected(1);
        end
    jPanelLeft.add(jPanelMeasureMN);
    
    % ==== PANEL: MEASURE DIPOLE MODELING ====
    jPanelMeasureDip = gui_river([1,1], [0,6,6,6], 'Measure');
        jGroupDipMeasure = ButtonGroup();
        jRadioDipBest    = gui_component('radio', jPanelMeasureDip, [],   'Best dipole fit',       jGroupDipMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioDipFit     = gui_component('radio', jPanelMeasureDip, 'br', 'Goodness-of-fit map',   jGroupDipMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioDipChi     = gui_component('radio', jPanelMeasureDip, 'br', 'Chi-squared error map', jGroupDipMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioDipPerform = gui_component('radio', jPanelMeasureDip, 'br', 'NP performance index',  jGroupDipMeasure, '', @(h,ev)UpdatePanel(), []);
        % Default selection
        switch lower(OPTIONS.InverseMeasure)
            case 'amplitude',    jRadioDipBest.setSelected(1);
            case 'fit',          jRadioDipFit.setSelected(1);
            case 'chi',          jRadioDipChi.setSelected(1);
            case 'performance',  jRadioDipPerform.setSelected(1);
            otherwise,           jRadioDipBest.setSelected(1);
        end
    jPanelLeft.add(jPanelMeasureDip);
    
    % ==== PANEL: MEASURE BEAMFORMER ====
    jPanelMeasureBf = gui_river([1,1], [0,6,6,6], 'Measure');
        jGroupBfMeasure = ButtonGroup();
        jRadioMethodBfTs      = gui_component('radio', jPanelMeasureBf, [],   'Beamformer time-series', jGroupBfMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioMethodBfPow     = gui_component('radio', jPanelMeasureBf, 'br', 'Beamformer power',       jGroupBfMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioMethodBfNai     = gui_component('radio', jPanelMeasureBf, 'br', 'Neural activity index',  jGroupBfMeasure, '', @(h,ev)UpdatePanel(), []);
        jRadioMethodBfPerform = gui_component('radio', jPanelMeasureBf, 'br', 'NP performance index',   jGroupBfMeasure, '', @(h,ev)UpdatePanel(), []);
        % Default selection
        switch lower(OPTIONS.InverseMeasure)
            case 'amplitude',    jRadioMethodBfTs.setSelected(1);
            case 'power',        jRadioMethodBfPow.setSelected(1);
            case 'nai',          jRadioMethodBfNai.setSelected(1);
            case 'performance',  jRadioMethodBfPerform.setSelected(1);
            otherwise,           jRadioMethodBfTs.setSelected(1);
        end
    jPanelLeft.add(jPanelMeasureBf);
    
    % ==== PANEL: SOURCE MODEL ====
    jPanelModel = gui_river([1,1], [0,6,6,6], 'Source model: Dipoles orientations');
        jGroupModel    = ButtonGroup(); 
        jRadioConstr   = gui_component('radio', jPanelModel, [],   'Constrained:  Normal to cortex',    jGroupModel, '', @(h,ev)UpdatePanel(), []);
        jRadioOptimal  = gui_component('radio', jPanelModel, 'br', 'Constrained:  Optimal orientation', jGroupModel, '', @(h,ev)UpdatePanel(), []);
        jRadioLoose    = gui_component('radio', jPanelModel, 'br', 'Loose constraints',                jGroupModel, '', @(h,ev)UpdatePanel(), []);
        jTextLoose     = gui_component('texttime', jPanelModel, [], '', [], '', [], []);
        gui_validate_text(jTextLoose, [], [], 0:0.1:1, '', 1, OPTIONS.Loose, []);
        jRadioUnconstr = gui_component('radio', jPanelModel, 'br', 'Unconstrained', jGroupModel, '', @(h,ev)UpdatePanel(), []);
        % Default selection
        if strcmpi(HeadModelType, 'surface')
            switch lower(OPTIONS.SourceOrient{1})
                case 'fixed',    jRadioConstr.setSelected(1);
                case 'optimal',  jRadioOptimal.setSelected(1);
                case 'loose',    jRadioLoose.setSelected(1);
                case 'free',     jRadioUnconstr.setSelected(1);
            end
        elseif strcmpi(HeadModelType, 'volume')
            jRadioConstr.setEnabled(0);
            jRadioLoose.setEnabled(0);
            jRadioUnconstr.setSelected(1);
        elseif strcmpi(HeadModelType, 'mixed')
            jRadioConstr.setEnabled(0);
            jRadioLoose.setEnabled(0);
            jRadioUnconstr.setEnabled(0);
        end
    jPanelLeft.add(jPanelModel);
    
    % ==== PANEL: DATA TYPE ====
    jPanelSensors = gui_river([1,1], [0,6,6,6], 'Sensors');
        jCheckMeg = [];
        jCheckMegGrad = [];
        jCheckMegMag = [];
        jCheckEeg = [];
        jCheckEcog = [];
        jCheckSeeg = [];
        if ismember('MEG', Modalities)
            jCheckMeg = gui_component('checkbox', jPanelSensors, '', 'MEG', [], '', @Modality_Callback, []);
            jCheckMeg.setSelected(1);
        end
        if ismember('MEG GRAD', Modalities)
            jCheckMegGrad = gui_component('checkbox', jPanelSensors, '', 'MEG GRAD', [], '', @Modality_Callback, []);
            jCheckMegGrad.setSelected(1);
        end
        if ismember('MEG MAG', Modalities)
            jCheckMegMag = gui_component('checkbox', jPanelSensors, '', 'MEG MAG', [], '', @Modality_Callback, []);
            jCheckMegMag.setSelected(1);
        end
        if ismember('EEG', Modalities)
            jCheckEeg = gui_component('checkbox', jPanelSensors, '', 'EEG', [], '', @Modality_Callback, []);
            if (length(Modalities) == 1)
                jCheckEeg.setSelected(1);
            end
        end
        if ismember('ECOG', Modalities)
            jCheckEcog = gui_component('checkbox', jPanelSensors, '', 'ECOG', [], '', @Modality_Callback, []);
            if (length(Modalities) == 1)
                jCheckEcog.setSelected(1);
            end
        end
        if ismember('SEEG', Modalities)
            jCheckSeeg = gui_component('checkbox', jPanelSensors, '', 'SEEG', [], '', @Modality_Callback, []);
            if (length(Modalities) == 1)
                jCheckSeeg.setSelected(1);
            end
        end
    jPanelLeft.add(jPanelSensors);
    
    % ======================================================================================================
    
    % ==== DEPTH WEIGHTING =====
    jPanelDepth = gui_river([1,1], [0,6,6,6], 'Depth weighting');
        % Use depth weighting
        jCheckDepth = gui_component('checkbox', jPanelDepth, 'br', 'Use depth weighting', [], '', @(h,ev)UpdatePanel(), []);
        jCheckDepth.setSelected(OPTIONS.UseDepth);
        % Weightexp
        jLabelWeightExp = gui_component('label', jPanelDepth, 'br', '       Order [0,1]: ', [], '', [], []);
        jTextWeightExp  = gui_component('texttime', jPanelDepth, 'tab', num2str(OPTIONS.WeightExp), [], '', [], []);
        % Weightlimit
        jLabelWeightLimit = gui_component('label', jPanelDepth, 'br', '       Maximal amount: ', [], '', [], []);
        jTextWeightLimit  = gui_component('texttime', jPanelDepth, 'tab', num2str(OPTIONS.WeightLimit), [], '', [], []);
    c.gridy = 1;
    jPanelRight.add(jPanelDepth, c);
    
    % ==== PANEL: NOISE COVARIANCE ====
    jPanelNoiseCov = gui_river([1,1], [0,6,6,6], 'Noise covariance regularization');
        jGroupReg    = ButtonGroup();
        jRadioShrink = gui_component('radio',    jPanelNoiseCov, [],   'Automatic shrinkage (Ledoit and Wolf)', jGroupReg, '', @(h,ev)UpdatePanel(), []);
        jRadioReg    = gui_component('radio',    jPanelNoiseCov, 'br', 'Regularize noise covariance', jGroupReg, '', @(h,ev)UpdatePanel(), []);
        jTextReg     = gui_component('texttime', jPanelNoiseCov, [], '', [], '', [], []);
        gui_validate_text(jTextReg, [], [], 0:0.1:1, '', 1, OPTIONS.NoiseReg, []);
        jRadioDiag   = gui_component('radio',    jPanelNoiseCov, 'br', 'Diagonal noise covariance', jGroupReg, '', @(h,ev)UpdatePanel(), []);
        jRadioNoReg  = gui_component('radio',    jPanelNoiseCov, 'br', 'None', jGroupReg, '', @(h,ev)UpdatePanel(), []);
        % Default selection
        switch lower(OPTIONS.NoiseMethod)
            case 'shrink', jRadioShrink.setSelected(1);
            case 'reg',    jRadioReg.setSelected(1);
            case 'diag',   jRadioDiag.setSelected(1);
            case 'none',   jRadioNoReg.setSelected(1);
        end
    c.gridy = 2;
    jPanelRight.add(jPanelNoiseCov, c);
    
    % ==== PANEL: SNR ====
    jPanelSnr = gui_river([1,1], [0,6,6,6], 'Signal-to-noise ratio');
        jGroupSnr  = ButtonGroup();
        % Maximum source amplitude
        jRadioSnrRms = gui_component('radio', jPanelSnr, [], 'RMS source amplitude: ', jGroupSnr, '', @(h,ev)UpdatePanel(), []);
        jTextSnrRms  = gui_component('texttime', jPanelSnr, [], '', [], '', [], []);
        gui_validate_text(jTextSnrRms, [], [], [0, 1000000, 100], 'scalar', 2, OPTIONS.SnrRms * 1e9, []);
        gui_component('label', jPanelSnr, [], 'nAm', [], '', [], []);
        % Fixed SNR
        jRadioSnrFix = gui_component('radio', jPanelSnr, 'br', 'Use fixed SNR: ', jGroupSnr, '', @(h,ev)UpdatePanel(), []);
        jTextSnrFix  = gui_component('texttime', jPanelSnr, [], '', [], '', [], []);
        gui_validate_text(jTextSnrFix, [], [], [0, 10000, 100], '', 2, OPTIONS.SnrFixed, []);
        jRadioSnrFix.setSelected(1);
        % Estimate SNR
        jRadioSnrEst = gui_component('radio', jPanelSnr, 'br', 'Estimate SNR from data', jGroupSnr, '', @(h,ev)UpdatePanel(), []);
        gui_component('label', jPanelSnr, 'br', '<HTML>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<FONT color="#555555"><I>SNR = max(eig(DataCov)) / max(eig(NoiseCov))</I></FONT>', [], '', [], []);
        % Default selection
        switch lower(OPTIONS.SnrMethod)
            case 'rms',      jRadioSnrRms.setSelected(1);
            case 'fixed',    jRadioSnrFix.setSelected(1);
            case 'estimate', jRadioSnrEst.setSelected(1);
        end
    c.gridy = 3;
    jPanelRight.add(jPanelSnr, c);
    
    % ===== PANEL: OUTPUT MODE =====
    jPanelOutput = gui_river([1,1], [0,6,6,6], 'Output mode');
        jGroupOutput = ButtonGroup();
        jRadioKernel = gui_component('radio', jPanelOutput, [],   'Inverse kernel only', jGroupOutput, '<HTML>Time independant computation.<BR>To get the sources estimations for a time frame, <BR> the kernel is applied to the recordings (matrix product).', @(h,ev)UpdatePanel(), []);
        jRadioFull   = gui_component('radio', jPanelOutput, 'br', 'Full results (Kernel*Recordings)', jGroupOutput, 'Compute sources for all the time samples.', @(h,ev)UpdatePanel(), []);
        if isShared
            jRadioFull.setEnabled(0);
            jRadioKernel.setSelected(1);
        elseif isFull
            jRadioFull.setSelected(1);
        else
            jRadioKernel.setSelected(1);
        end
    c.gridy = 4;
    jPanelRight.add(jPanelOutput, c);
    
    % ===== VALIDATION BUTTONS =====
    jPanelValid = gui_river([1,1], [0,6,6,6]);
    % Expert/normal mode
    jButtonExpert = gui_component('button', jPanelValid, [], 'Show details', [], [], @SwitchExpertMode_Callback, []);
    gui_component('label', jPanelValid, 'hfill', ' ');
    % Ok/Cancel
    gui_component('Button', jPanelValid, 'right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('Button', jPanelValid, [], 'OK', [], [], @ButtonOk_Callback, []);
    jPanelLeft.add(jPanelValid);
    % Add a glue at the bottom of the right panel (for appropriate scaling with the left)
    c.gridy   = 5;
    c.weighty = 1;
    jPanelRight.add(Box.createVerticalGlue(), c);

    % ===== DISABLE EXPERIMENTAL WORK =====
    % Hide process if not explicitly allow
    if isempty(ALLOW_2015) || ~ALLOW_2015
        jRadioMnSloreta.setEnabled(0);
        jRadioDipFit.setEnabled(0);
        jRadioDipChi.setEnabled(0);
        jRadioMethodBfPow.setEnabled(0);
        jRadioMethodBfNai.setEnabled(0);
        jRadioOptimal.setEnabled(0);
        jRadioSnrEst.setEnabled(0);
        jRadioShrink.setEnabled(0);
        if jRadioShrink.isSelected()
            jRadioReg.setSelected(1);
        end
    end
 
    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Create the BstPanel object that is returned by the function
    ctrl = struct(...
            'HeadModelType',  HeadModelType, ...
            'jTextComment',   jTextComment, ...
            'jToogleLinear',  jToogleLinear, ...
            'jToogleNonLin',  jToogleNonLin, ...
            ... % ==== PANEL: METHOD ====
            'jRadioMethodMn',  jRadioMethodMn, ...
            'jRadioMethodBf',  jRadioMethodBf, ...
            'jRadioMethodDip', jRadioMethodDip, ...
            ... % ==== PANEL: MEASURE ====
            'jRadioMnCurrent', jRadioMnCurrent, ...
            'jRadioMnDspm',    jRadioMnDspm, ...
            'jRadioMnPerform', jRadioMnPerform, ...
            'jRadioMnSloreta', jRadioMnSloreta, ...
            'jRadioDipBest',     jRadioDipBest, ...
            'jRadioDipFit',      jRadioDipFit, ...
            'jRadioDipChi',      jRadioDipChi, ...
            'jRadioDipPerform',  jRadioDipPerform, ...
            'jRadioMethodBfTs',      jRadioMethodBfTs, ...
            'jRadioMethodBfPow',     jRadioMethodBfPow, ...
            'jRadioMethodBfNai',     jRadioMethodBfNai, ...
            'jRadioMethodBfPerform', jRadioMethodBfPerform, ...
            ... % ==== PANEL: SOURCE MODEL ====
            'jRadioConstr',   jRadioConstr, ...
            'jRadioOptimal',  jRadioOptimal, ...
            'jRadioUnconstr', jRadioUnconstr, ...
            'jRadioLoose',    jRadioLoose, ...
            'jTextLoose',     jTextLoose, ...
            ... % ==== PANEL: DEPTH WEIGHTING ====
            'jCheckDepth',      jCheckDepth, ...
            'jTextWeightExp',   jTextWeightExp, ...
            'jTextWeightLimit', jTextWeightLimit, ...
            ... % ==== PANEL: SNR ====
            'jRadioSnrRms',  jRadioSnrRms, ...
            'jTextSnrRms',   jTextSnrRms, ...
            'jRadioSnrFix',  jRadioSnrFix, ...
            'jTextSnrFix',   jTextSnrFix, ...
            'jRadioSnrEst',  jRadioSnrEst, ...
            ... % ==== PANEL: NOISE COVARIANCE ====
            'jRadioShrink',  jRadioShrink, ...
            'jRadioReg',     jRadioReg, ...
            'jTextReg',      jTextReg, ...
            'jRadioDiag',    jRadioDiag, ...
            'jRadioNoReg',   jRadioNoReg, ...
            ... % ==== PANEL: NON-LINEAR ====
            'jRadioMem',     jRadioMem, ...
            ... % ==== PANEL: DATA TYPE ====
            'jCheckMeg',     jCheckMeg, ...
            'jCheckMegGrad', jCheckMegGrad, ...
            'jCheckMegMag',  jCheckMegMag, ...
            'jCheckEeg',     jCheckEeg, ...
            'jCheckEcog',    jCheckEcog, ...
            'jCheckSeeg',    jCheckSeeg, ...
            ... % ===== PANEL: OUTPUT MODE =====
            'jRadioFull',    jRadioFull, ...
            'jRadioKernel',  jRadioKernel);
    % Create the BstPanel object that is returned by the function
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
    % Update comments
    UpdatePanel(1, 1);
    


%% =================================================================================
%  === LOCAL CALLBACKS  ============================================================
%  =================================================================================
    %% ===== BUTTON: CANCEL =====
    function ButtonCancel_Callback(varargin)
        % Close panel
        gui_hide(panelName);
    end

    %% ===== BUTTON: OK =====
    function ButtonOk_Callback(varargin)
        % Release mutex and keep the panel opened
        bst_mutex('release', panelName);
    end

    %% ===== MODALITY CALLBACK =====
    function Modality_Callback(hObject, event)
        % If only one checkbox: can't deselect it
        if (length(Modalities) == 1)
            event.getSource().setSelected(1);
        % Warning if both MEG and EEG are selected
        elseif isFirstCombinationWarning && ~isempty(jCheckEeg) && jCheckEeg.isSelected() && (...
                (~isempty(jCheckMeg) && jCheckMeg.isSelected()) || ...
                (~isempty(jCheckMegGrad) && jCheckMegGrad.isSelected()) || ...
                (~isempty(jCheckMegMag) && jCheckMegMag.isSelected()))
            java_dialog('warning', ['Warning: Brainstorm inverse models do not properly handle the combination of MEG and EEG yet.' 10 10 ...
                                    'For now, we recommend to compute separatly the sources for MEG and EEG.'], 'EEG/MEG combination');
            isFirstCombinationWarning = 0;
        end
        % Update comment
        UpdatePanel();
    end


    %% ===== SWITCH EXPERT MODE =====
    function SwitchExpertMode_Callback(varargin)
        % Toggle expert mode
        ExpertMode = bst_get('ExpertMode');
        bst_set('ExpertMode', ~ExpertMode);
        % Update comment
        UpdatePanel(1);
    end
    

    %% ===== UPDATE PANEL ======
    % USAGE:  UpdatePanel(isForced = 0)
    function UpdatePanel(isForced, isFirstCall)
        % Default values
        if (nargin < 2) || isempty(isFirstCall)
            isFirstCall = 0;
        end
        if (nargin < 1) || isempty(isForced)
            isForced = 0;
        end
        % Get the main categories of options
        isLinear = jToogleLinear.isSelected();
        % Expert mode / Normal mode
        if isForced
            ExpertMode = bst_get('ExpertMode');
            % Show/hide left panels
            jPanelMethod.setVisible(isLinear);
            % jPanelModel.setVisible(isLinear && ~strcmpi(HeadModelType, 'mixed'));
            jPanelModel.setVisible(isLinear);
            jPanelMeasureMN.setVisible(isLinear  && jRadioMethodMn.isSelected());
            jPanelMeasureDip.setVisible(isLinear && jRadioMethodDip.isSelected());
            jPanelMeasureBf.setVisible(isLinear  && jRadioMethodBf.isSelected());
            jPanelNonLin.setVisible(~isLinear);
            % Right panels (expert)
            jPanelRight.setVisible(ExpertMode);
            jPanelNoiseCov.setVisible(isLinear);
            % jPanelSnr.setVisible(isLinear && jRadioMethodMn.isSelected());
            jPanelSnr.setVisible(isLinear);
            jPanelDepth.setVisible(isLinear);
            jPanelOutput.setVisible(isLinear && ~isProcess);
            % Update expert button 
            if ExpertMode
                jButtonExpert.setText('Hide details');
            else
                jButtonExpert.setText('Show details');
            end
            % Get old panel
            [bstPanelOld, iPanel] = bst_get('Panel', 'InverseOptions');
            container = get(bstPanelOld, 'container');
            % Re-pack frame
            if ~isempty(container)
                jFrame = container.handle{1};
                if ~isempty(jFrame)
                    jFrame.pack();
                end
            end
        end
        % Enable/disable text boxes
        jTextLoose.setEnabled(jRadioLoose.isSelected());
        jLabelWeightExp.setEnabled(jCheckDepth.isSelected());
        jTextWeightExp.setEnabled(jCheckDepth.isSelected());
        jLabelWeightLimit.setEnabled(jCheckDepth.isSelected());
        jTextWeightLimit.setEnabled(jCheckDepth.isSelected());
        jTextSnrRms.setEnabled(jRadioSnrRms.isSelected());
        jTextSnrFix.setEnabled(jRadioSnrFix.isSelected());
        jTextReg.setEnabled(jRadioReg.isSelected());
        % Selected modalities
        selModalities = {};
        if ~isempty(jCheckMeg) && jCheckMeg.isSelected()
            selModalities{end+1} = 'MEG';
        end
        if ~isempty(jCheckMegGrad) && jCheckMegGrad.isSelected()
            selModalities{end+1} = 'MEG GRAD';
        end
        if ~isempty(jCheckMegMag) && jCheckMegMag.isSelected()
            selModalities{end+1} = 'MEG MAG';
        end
        if ~isempty(jCheckEeg) && jCheckEeg.isSelected()
            selModalities{end+1} = 'EEG';
        end
        if ~isempty(jCheckEcog) && jCheckEcog.isSelected()
            selModalities{end+1} = 'ECOG';
        end
        if ~isempty(jCheckSeeg) && jCheckSeeg.isSelected()
            selModalities{end+1} = 'SEEG';
        end
        
        % Linear methods
        if isLinear
            % Get selected method
            [InverseMethod, InverseMeasure] = GetSelectedMethod(ctrl);
            % Get comment for this method
            Comment = GetMethodComment(InverseMethod, InverseMeasure);
        else
            if jRadioMem.isSelected()
                Comment = 'MEM: ';
            end
        end
        % Add modality comment
        Comment = [Comment, ': ', process_inverse_2015('GetModalityComment', selModalities)];
        % Update comment field
        if ~isProcess || ~isFirstCall || ~isfield(OPTIONS, 'Comment') || isempty(OPTIONS.Comment)
            jTextComment.setText(Comment);
        end
    end
end


%% =================================================================================
%  === EXTERNAL CALLBACKS  =========================================================
%  =================================================================================
%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    % Get panel controls handles
    ctrl = bst_get('PanelControls', 'InverseOptions');
    if isempty(ctrl)
        s = [];
        return; 
    end
    % Comment
    s.Comment = char(ctrl.jTextComment.getText());
    % Linear models
    if ctrl.jToogleLinear.isSelected()
        % Get selected method
        [s.InverseMethod, s.InverseMeasure] = GetSelectedMethod(ctrl);
        % Source model
        if strcmpi(ctrl.HeadModelType, 'mixed')
            s.SourceOrient = [];
            s.Loose = [];
        else
            if ctrl.jRadioConstr.isSelected()
                s.SourceOrient = {'fixed'};
            elseif ctrl.jRadioOptimal.isSelected()
                s.SourceOrient = {'optimal'};
            elseif ctrl.jRadioUnconstr.isSelected()
                s.SourceOrient = {'free'};
            elseif ctrl.jRadioLoose.isSelected()
                s.SourceOrient = {'loose'};
            end
            s.Loose = str2num(char(ctrl.jTextLoose.getText()));
        end
        % Depth weighting
        s.UseDepth    = ctrl.jCheckDepth.isSelected();
        s.WeightExp   = str2num(char(ctrl.jTextWeightExp.getText()));
        s.WeightLimit = str2num(char(ctrl.jTextWeightLimit.getText()));
        % Noise covariance
        if ctrl.jRadioShrink.isSelected()
            s.NoiseMethod = 'shrink';
        elseif ctrl.jRadioReg.isSelected()
            s.NoiseMethod = 'reg';
        elseif ctrl.jRadioDiag.isSelected()
            s.NoiseMethod = 'diag';
        elseif ctrl.jRadioNoReg.isSelected()
            s.NoiseMethod = 'none';
        end
        s.NoiseReg = str2num(char(ctrl.jTextReg.getText()));
        % Signal to noise
        if ctrl.jRadioSnrRms.isSelected()
            s.SnrMethod = 'rms';
        elseif ctrl.jRadioSnrFix.isSelected()
            s.SnrMethod = 'fixed';
        elseif ctrl.jRadioSnrEst.isSelected()
            s.SnrMethod = 'estimate';
        end
        s.SnrRms   = str2num(char(ctrl.jTextSnrRms.getText())) * 1e-9;  % Convert to Amper.mter
        s.SnrFixed = str2num(char(ctrl.jTextSnrFix.getText()));
        % Output mode
        if ctrl.jRadioFull.isSelected()
            s.ComputeKernel = 0;
        else
            s.ComputeKernel = 1;
        end
    % Non-linear models
    else
        % Get selected method
        if ctrl.jRadioMem.isSelected()
            s.InverseMethod = 'mem';
        end
        % Other fields that are not defined
        s.InverseMeasure = [];
        s.SourceOrient   = [];
        s.Loose          = [];
        s.UseDepth       = [];
        s.WeightExp      = [];
        s.WeightLimit    = [];
        s.NoiseMethod    = [];
        s.NoiseReg       = [];
        s.SnrMethod      = [];
        s.SnrRms         = [];
        s.SnrFixed       = [];
        s.ComputeKernel  = 0;
    end
    % Selected modalities
    s.DataTypes = {};
    if ~isempty(ctrl.jCheckMeg) && ctrl.jCheckMeg.isSelected()
        s.DataTypes{end+1} = 'MEG';
    end
    if ~isempty(ctrl.jCheckMegGrad) && ctrl.jCheckMegGrad.isSelected()
        s.DataTypes{end+1} = 'MEG GRAD';
    end
    if ~isempty(ctrl.jCheckMegMag) && ctrl.jCheckMegMag.isSelected()
        s.DataTypes{end+1} = 'MEG MAG';
    end
    if ~isempty(ctrl.jCheckEeg) && ctrl.jCheckEeg.isSelected()
        s.DataTypes{end+1} = 'EEG';
    end
    if ~isempty(ctrl.jCheckEcog) && ctrl.jCheckEcog.isSelected()
        s.DataTypes{end+1} = 'ECOG';
    end
    if ~isempty(ctrl.jCheckSeeg) && ctrl.jCheckSeeg.isSelected()
        s.DataTypes{end+1} = 'SEEG';
    end
end


%% ===== GET SELECTED METHOD =====
function [Method, Measure] = GetSelectedMethod(ctrl)
    if ctrl.jRadioMethodMn.isSelected()
        Method = 'minnorm';
        if ctrl.jRadioMnCurrent.isSelected()
            Measure = 'amplitude';
        elseif ctrl.jRadioMnPerform.isSelected()
            Measure = 'performance';
        elseif ctrl.jRadioMnDspm.isSelected()
            Measure = 'dspm';
        elseif ctrl.jRadioMnSloreta.isSelected()
            Measure = 'sloreta';
        end
    elseif ctrl.jRadioMethodDip.isSelected()
        Method = 'gls';
        if ctrl.jRadioDipBest.isSelected()
            Measure = 'amplitude';
        elseif ctrl.jRadioDipFit.isSelected()
            Measure = 'fit';
        elseif ctrl.jRadioDipChi.isSelected()
            Measure = 'chi';
        elseif ctrl.jRadioDipPerform.isSelected()
            Measure = 'performance';
        end
    elseif ctrl.jRadioMethodBf.isSelected()
        Method = 'lcmv';
        if ctrl.jRadioMethodBfTs.isSelected()
            Measure = 'amplitude';
        elseif ctrl.jRadioMethodBfPow.isSelected()
            Measure = 'power';
        elseif ctrl.jRadioMethodBfNai.isSelected()
            Measure = 'nai';
        elseif ctrl.jRadioMethodBfPerform.isSelected()
            Measure = 'performance';
        end
    end
end

%% ===== GET METHOD COMMENT =====
function Comment = GetMethodComment(Method, Measure)
    switch (lower(Method))
        case 'minnorm'
            switch (lower(Measure))
                case 'amplitude',    Comment = 'MN';
                case 'performance',  Comment = 'MNp';
                case 'dspm',         Comment = 'dSPM';
                case 'sloreta',      Comment = 'sLORETA';
            end
        case 'gls'
            switch (lower(Measure))
                case 'amplitude',    Comment = 'GLS';
                case 'fit',          Comment = 'GLSfit';
                case 'chi',          Comment = 'GLSchi';
                case 'performance',  Comment = 'GLSp';
            end
        case 'lcmv'
            switch (lower(Measure))
                case 'amplitude',    Comment = 'LCMV';
                case 'power',        Comment = 'LCMVpow';
                case 'nai',          Comment = 'NAI';
                case 'performance',  Comment = 'LCMVp';
            end
        case 'mem'
            Comment = 'MEM';
    end
end

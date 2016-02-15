function varargout = process_bandpass( varargin )
% PROCESS_BANDPASS: Frequency filters: Lowpass/Highpass/Bandpass
%
% USAGE:      sProcess = process_bandpass('GetDescription')
%               sInput = process_bandpass('Run', sProcess, sInput, method=[])
%                    x = process_bandpass('Compute', x, sfreq, HighPass, LowPass, method=[], isMirror=1)

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
% Authors: Francois Tadel, 2010-2014

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Band-pass filter';
    sProcess.FileTag     = '| bandpass';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 64;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'raw', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'raw', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.processDim  = 1;   % Process channel by channel
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/Epilepsy#Band-pass_filter';
    % Definition of the options
    % === Low bound
    sProcess.options.highpass.Comment = 'Lower cutoff frequency (0=disable):';
    sProcess.options.highpass.Type    = 'value';
    sProcess.options.highpass.Value   = {2,'Hz ',2};
    % === High bound
    sProcess.options.lowpass.Comment = 'Upper cutoff frequency (0=disable):';
    sProcess.options.lowpass.Type    = 'value';
    sProcess.options.lowpass.Value   = {40,'Hz ',2};
    % === Mirror
    sProcess.options.mirror.Comment = 'Mirror signal before filtering (to avoid edge effects)';
    sProcess.options.mirror.Type    = 'checkbox';
    sProcess.options.mirror.Value   = 1;
    % === Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data', 'raw'};
end

%% ===== GET OPTIONS =====
function [HighPass, LowPass, isMirror] = GetOptions(sProcess)
    HighPass = sProcess.options.highpass.Value{1};
    LowPass  = sProcess.options.lowpass.Value{1};
    if (HighPass == 0) 
        HighPass = [];
    end
    if (LowPass == 0) 
        LowPass = [];
    end
    isMirror = sProcess.options.mirror.Value;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    % Get options
    [HighPass, LowPass] = GetOptions(sProcess);
    % Format comment
    if ~isempty(HighPass) && ~isempty(LowPass)
        Comment = ['Band-pass:' num2str(HighPass) 'Hz-' num2str(LowPass) 'Hz'];
    elseif ~isempty(HighPass)
        Comment = ['High-pass:' num2str(HighPass) 'Hz'];
    elseif ~isempty(LowPass)
        Comment = ['Low-pass:' num2str(LowPass) 'Hz'];
    else
        Comment = '';
    end
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
    % Get options
    [HighPass, LowPass, isMirror] = GetOptions(sProcess);
    % ===== FILTER DATA =====
    sfreq = 1 ./ (sInput.TimeVector(2) - sInput.TimeVector(1));
    sInput.A = Compute(sInput.A, sfreq, HighPass, LowPass, [], isMirror);
    % ===== FILE COMMENT =====
    if ~isempty(HighPass) && ~isempty(LowPass)
        filterComment = ['| band(' num2str(HighPass) '-' num2str(LowPass) 'Hz)'];
    elseif ~isempty(HighPass)
        filterComment = ['| high(' num2str(HighPass) 'Hz)'];
    elseif ~isempty(LowPass)
        filterComment = ['| low(', num2str(LowPass) 'Hz)'];
    else
        filterComment = '';
    end
    sInput.FileTag = filterComment;
    % Comment for History
    sInput.HistoryComment = strrep(filterComment, '| ', '');
    % Do not keep the Std field in the output
    if isfield(sInput, 'Std') && ~isempty(sInput.Std)
        sInput.Std = [];
    end
end


%% ===== EXTERNAL CALL =====
% USAGE: process_bandpass('Compute', x, sfreq, HighPass, LowPass, method=[], isMirror=1)
function x = Compute(x, sfreq, HighPass, LowPass, method, isMirror)
    % Default method
    if (nargin < 6) || isempty(isMirror)
        isMirror = 1;
    end
    if (nargin < 5) || isempty(method)
        method = 'bst-fft-fir';
    end
    % Filtering using the selected method
    switch (method)
        % Better filter, a bit slower
        % WARNING: APPLIES A LOW-PASS FILTER Fs/3 IF WE WANT A HIGH-PASS ONLY
        case 'bst-fft-fir'
            x = bst_bandpass_fft(x, sfreq, HighPass, LowPass, 1, isMirror);
        % Doesn't make a big difference, faster, but not all the cases are handled properly
        % WARNING: APPLIES A LOW-PASS FILTER Fs/3 IF WE WANT A HIGH-PASS ONLY
        case 'bst-fft'
            x = bst_bandpass_fft(x, sfreq, HighPass, LowPass, 0, isMirror);
        % Filter using SOS functions: too slow, unstable...
        case 'bst-sos'
            % Prepare options structure
            coef.LowPass = LowPass;
            coef.HighPass = HighPass;
            % Filter signal
            x = bst_bandpass_sos(x, sfreq, coef);
            
        % New filters designed by JC Mosher
        case 'bst-filtfilt-fir'
            x = bst_bandpass_filtfilt(x, sfreq, HighPass, LowPass, 0, 'fir');
        case 'bst-filtfilt-iir'
            x = bst_bandpass_filtfilt(x, sfreq, HighPass, LowPass, 0, 'iir');
        
        % Source: FieldTrip toolbox
        % Equivalent to: x = ft_preproc_bandpassfilter(x, sfreq, FreqBand, [], 'but');
        case 'fieldtrip_butter'
            % Nyqist frequency
            Fnyq = sfreq/2;
            FreqBand = [HighPass, LowPass];
            % Filter order
            N = 4;
            % Butterworth filter
            if bst_get('UseSigProcToolbox')
                [B,A] = butter(N, FreqBand ./ Fnyq);
            else
                [B,A] = oc_butter(N, FreqBand ./ Fnyq);
            end
            % Filter
            x = filtfilt(B, A, x')';
            

        % Source: FieldTrip toolbox
        % Bandpass filter: Onepass-zerophase, hamming-windowed sinc FIR
        % Equivalent to: x = ft_preproc_bandpassfilter(x, sfreq, FreqBand, [], 'firws');
        case 'fieldtrip_firws'
            % Nyqist frequency
            Fnyq = sfreq/2;
            FreqBand = [HighPass, LowPass];
            % Constants
            TRANSWIDTHRATIO = 0.25;
            % Max possible transition band width
            maxTBWArray = [FreqBand * 2, (Fnyq - FreqBand) * 2, diff(FreqBand)];
            maxDf = min(maxTBWArray);
            % Default filter order heuristic
            df = min([max([FreqBand(1) * TRANSWIDTHRATIO, 2]) maxDf]);
            if (df > maxDf)
                error('Transition band too wide. Maximum transition width is %.2f Hz.', maxDf)
            end
            % Compute filter order from transition width
            N = firwsord('hamming', sfreq, df, []);
            % Window
            win = bst_window('hamming', N+1);
            % Impulse response
            B = firws(N, FreqBand / Fnyq, 'band', win);
            % Padding
            x = x';
            groupDelay = (length(B) - 1) / 2;
            startPad = repmat(x(1,:), [groupDelay 1]);
            endPad = repmat(x(end,:), [groupDelay 1]);
            % Filter data
            x = filter(B, 1, [startPad; x; endPad]);
            % Remove padded data
            x = x(2 * groupDelay + 1:end, :);
            x = x';
    end
end




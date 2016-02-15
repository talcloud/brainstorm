function varargout = process_test_parametric2( varargin )
% PROCESS_TEST_PARAMETRIC2: Parametric independent two-sample tests.
% 
% USAGE:  OutputFiles = process_test_parametric2('Run', sProcess, sInput)
%                   p = process_test_parametric2('ComputePvalues', t, df, TestType='t', TestTail='two')

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
% Authors: Francois Tadel, Dimitrios Pantazis, 2008-2016
macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Parametric test: Independent';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Stat2';
    sProcess.SubGroup    = 'Test';
    sProcess.Index       = 101;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/Statistics';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data',  'results',  'timefreq',  'matrix'};
    sProcess.OutputTypes = {'pdata', 'presults', 'ptimefreq', 'pmatrix'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 2;

    % === GENERIC EXTRACT OPTIONS
    % Label
    sProcess.options.extract_title.Comment    = '<B><U>Select data to test</U></B>:';
    sProcess.options.extract_title.Type       = 'label';
    % Options
    sProcess = process_extract_values('DefineExtractOptions', sProcess);
    % Rename some of the default options
    sProcess.options.isnorm.Comment = 'Test absolute values (or norm for unconstrained sources)';
    sProcess.options.isabs.Comment  = 'Test absolute values';
    
    % === MATCH ROWS WITH NAMES
    sProcess.options.matchrows.Comment    = 'Match signals between files using their names';
    sProcess.options.matchrows.Type       = 'checkbox';
    sProcess.options.matchrows.Value      = 1;
    sProcess.options.matchrows.InputTypes = {'timefreq', 'matrix'};
    % === OUTPUT COMMENT
    sProcess.options.Comment.Comment = 'Comment (empty=default): ';
    sProcess.options.Comment.Type    = 'text';
    sProcess.options.Comment.Value   = '';
    
    % === TEST: title
    sProcess.options.test_title.Comment    = '<BR><B><U>Test function</U></B>:';
    sProcess.options.test_title.Type       = 'label';
    % === TEST: type
    sProcess.options.test_type.Comment = {['<B>Student''s t-test &nbsp;&nbsp;(equal variance)</B> <BR>t = (mean(A)-mean(B)) / (Sx * sqrt(1/nA + 1/nB))<BR>' ...
                                           'Sx = sqrt(((nA-1)*var(A) + (nB-1)*var(B)) / (nA+nB-2)) <BR>' ...
                                           '<FONT COLOR="#777777">df = nA + nB - 2</FONT>'], ...
                                          ['<B>Student''s t-test &nbsp;&nbsp;(unequal variance)</B> <BR>', ...
                                           't = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)<BR>' ...
                                           '<FONT COLOR="#777777">df=(vA/nA+vB/nB)<SUP>2</SUP> / ((vA/nA)<SUP>2</SUP>/(nA-1)+(vB/nB)<SUP>2</SUP>/(nB-1))</FONT>']; ...
                                          'ttest_equal', 'ttest_unequal'};
%     sProcess.options.test_type.Comment = {['<B>Student''s t-test &nbsp;&nbsp;(equal variance)</B> <BR>t = (mean(A)-mean(B)) / (Sx * sqrt(1/nA + 1/nB))<BR>' ...
%                                            'Sx = sqrt(((nA-1)*var(A) + (nB-1)*var(B)) ./ (nA+nB-2))'], ...
%                                           ['<B>Student''s t-test &nbsp;&nbsp;(unequal variance)</B>:<BR>', ...
%                                            't = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)'], ...
%                                           ['<B>Absolute mean test (zero-mean, half-normal)</B><BR>', ...
%                                            '???'], ...
%                                           ['<B>Absolute mean test (non zero-mean, folded-normal)</B>:<BR>', ...
%                                            '???'], ...
%                                           ['<B>Power test</B>:<BR>', ...
%                                            'F = (sum(Ai^2)/nA) / (sum(Bi^2)/nB) &nbsp;&nbsp;&nbsp;&nbsp; <I>~F(nA,nB)</I>'], ...
%                                           ['<B>Power test (unconstrained sources)</B>:<BR>', ...
%                                            '???']; ...
%                                           'ttest_equal', 'ttest_unequal', 'abstest', 'abstest2', 'power', 'power_unconstr'};
    sProcess.options.test_type.Type    = 'radio_label';
    sProcess.options.test_type.Value   = 'ttest_equal';
    % === TAIL FOR THE TEST STATISTIC
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'two';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    % === DATA SELECTION COMMENT ===
    strData = process_extract_time('GetTimeString', sProcess);
    if isfield(sProcess.options, 'freqrange') && isfield(sProcess.options.freqrange, 'Value') && iscell(sProcess.options.freqrange.Value) && (length(sProcess.options.freqrange.Value) == 3) && (length(sProcess.options.freqrange.Value{1}) == 2)
        FreqRange = sProcess.options.freqrange.Value{1};
        if (FreqRange(1) == FreqRange(2))
            strData = [strData, ' ' num2str(FreqRange(1)) 'Hz'];
        else
            strData = [strData, ' ' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) 'Hz'];
        end
    end
    if isfield(sProcess.options, 'sensortypes') && isfield(sProcess.options.sensortypes, 'Value') && ~isempty(sProcess.options.sensortypes.Value)
        strData = [strData, ' ', sProcess.options.sensortypes.Value];
    end
    if isfield(sProcess.options, 'rows') && isfield(sProcess.options.rows, 'Value') && ~isempty(sProcess.options.rows.Value)
        strData = [strData, ' ', sProcess.options.rows.Value];
    end
    if ~isempty(strData)
        strData = [' [' strData ']'];
    end

    % === ABSOLUTE VALUE ===
    % Get options
    if isfield(sProcess.options, 'isabs') && isfield(sProcess.options.isabs, 'Value') && sProcess.options.isabs.Value
        isAbsolute = 1;
        strAbs = ' abs';
    elseif isfield(sProcess.options, 'isnorm') && isfield(sProcess.options.isnorm, 'Value') && sProcess.options.isnorm.Value
        isAbsolute = 1;
        strAbs = ' norm';

    else
        isAbsolute = 0;
        strAbs = '';
    end
    
    % === TEST COMMENT ===
    % Get test info
    TestType = sProcess.options.test_type.Value;
    if isfield(sProcess.options, 'tail') && isfield(sProcess.options.tail, 'Value') && ~isempty(sProcess.options.tail.Value)
        TestTail = sProcess.options.tail.Value;
    else
        TestTail = [];
    end
    % Documenting test to perform
    switch (TestType)
        case {'ttest_equal', 'ttest_unequal', 'ttest_paired'}
            if isAbsolute
                strHypo = '          H0:( |A| - |B| = 0 )';
            else
                strHypo = '          H0:(A-B = 0)';
            end
        case 'ttest_abspaired'
            if isAbsolute
                strHypo = '          H0:( ||A| - |B|| = 0 )';
            else
                strHypo = '          H0:( |A-B| = 0 )';
            end
        case 'ttest_onesample'
            if isAbsolute
                strHypo = '          H0:( |X| = 0 )';
            else
                strHypo = '          H0:(X = 0)';
            end
        case 'chi2_onesample'
            strHypo = '          H0:(|Zi| = 0)';
        case {'abstest', 'abstest2'}
            strHypo = '          H0:(|mean(A)|-|mean(B)| = 0)';
    end
    strStat = 'test';
    
    % No comment when forcing a one-sided test 
    if ismember(TestType, {'ttest_abspaired', 'ttest_onesample'}) && isAbsolute && ismember(TestTail, {'one-', 'two', 'one+'})
        strTail = '';
    elseif strcmpi(TestType, 'chi2_onesample')
        strTail = '';
        strAbs = '';
    % Comment for no-test option (just computing the statistic)
    elseif isequal(TestTail, 'no')
        strStat = 'stat';
        strTail = '';
        strData = '';
        strHypo = '';
    % Comment for one-tailed tests
    elseif ismember(TestTail, {'one-','one+'})
        strTail = [' ' TestTail];
    else
        strTail = '';
    end

    % === ASSEMBLING ===
    switch (TestType)
        case 'ttest_equal',      Comment = ['t-' strStat ' [equal'   strTail strAbs ']' strData strHypo];
        case 'ttest_unequal',    Comment = ['t-' strStat ' [unequal' strTail strAbs ']' strData strHypo];
        case 'ttest_paired',     Comment = ['t-' strStat ' [paired'  strTail strAbs ']' strData strHypo];
        case 'ttest_abspaired',  Comment = ['t-' strStat ' [abspaired' strTail strAbs ']' strData strHypo];
        case 'ttest_onesample',  Comment = ['t-' strStat ' [zero' strTail strAbs ']' strData strHypo];
        case 'chi2_onesample',   Comment = ['Chi2-' strStat ' [zero' strTail ']' strData strHypo];
        case 'abstest',          Comment = ['Absolute test' strData strHypo];
        case 'abstest2',         Comment = ['Absolute test' strData strHypo];
        case 'power',            Comment = ['Power test' strData strHypo];
        case 'power_unconstr',   Comment = ['Power test (unconstrained)' strData strHypo];
    end
end


%% ===== RUN =====
function sOutput = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>
    % Initialize returned variables
    sOutput = [];
    
    % ===== GET OPTIONS =====
    % Get generic extract options
    OPTIONS = process_extract_values('GetExtractOptions', sProcess, sInputsA(1));
    % Match signals between files using their names
    if isfield(sProcess.options, 'matchrows') && isfield(sProcess.options.matchrows, 'Value') && ~isempty(sProcess.options.matchrows.Value)
        OPTIONS.isMatchRows = sProcess.options.matchrows.Value;
    else
        OPTIONS.isMatchRows = 1;
    end
    % Get test type
    TestType = sProcess.options.test_type.Value;
    TestTail = sProcess.options.tail.Value;
    % Invalid test/tail combinations
    if ismember(TestType, {'ttest_abspaired', 'ttest_onesample'}) && OPTIONS.isAbsolute && ismember(TestTail, {'two', 'one-'})
        bst_report('Warning', sProcess, [], 'Testing |X|>0: Using a positive one-tailed test (one+) instead.');
        TestTail = 'one+';
    elseif strcmpi(TestType, 'chi2_onesample') && ismember(TestTail, {'two', 'one-'})
        bst_report('Warning', sProcess, [], 'Testing |X|>0: Using a positive one-tailed test (one+) instead.');
        TestTail = 'one+';
    end
    % Get average function
    switch (TestType)
        case {'ttest_equal', 'ttest_unequal', 'ttest_onesample', 'ttest_paired', 'abstest', 'abstest2'}
            if OPTIONS.isAbsolute
                AvgFunction = 'norm';
                isAvgVariance = 1;
            else
                AvgFunction = 'mean';
                isAvgVariance = 1;
            end
        case 'ttest_abspaired'
            if OPTIONS.isAbsolute
                AvgFunction = 'normdiffnorm';
            else
                AvgFunction = 'meandiffnorm';
            end
            isAvgVariance = 1;
        case {'power', 'power_unconstr', 'chi2_onesample'}
            AvgFunction = 'rms';
            isAvgVariance = 0;
    end
    isAvgWeighted = 0;
            
    % ===== CHECK INPUT FILES =====
    % Make sure that file type is indentical for both sets
    if ~isempty(sInputsA) && ~isempty(sInputsB) && ~strcmpi(sInputsA(1).FileType, sInputsB(1).FileType)
        bst_report('Error', sProcess, [], 'Cannot process inputs from different types.');
        return;
    end
    % Check the number of files in input
    if (length(sInputsA) < 2)
        bst_report('Error', sProcess, [], 'Not enough files in input.');
        return;
    end
    % Load time vector from the first file: if same as input, discard input
    TimeVector = in_bst(sInputsA(1).FileName, 'Time');
    if ~isempty(OPTIONS.TimeWindow) && (abs(TimeVector(1) - OPTIONS.TimeWindow(1)) < 1e-4) && (abs(TimeVector(end) - OPTIONS.TimeWindow(2)) < 1e-4)
        OPTIONS.TimeWindow = [];
    end
    % Load freq range from the first file: if same as input, discard input
    if ~isempty(OPTIONS.FreqRange)
        % Load Freqs field from the input file
        TfMat = in_bst_timefreq(sInputsA(1).FileName, 0, 'Freqs');
        if iscell(TfMat.Freqs)
            BandBounds = process_tf_bands('GetBounds', TfMat.Freqs);
            FreqList = unique(BandBounds(:));
        else
            FreqList = TfMat.Freqs;
        end
        if (OPTIONS.FreqRange(1) - FreqList(1) < 1e-4) && (OPTIONS.FreqRange(2) - FreqList(end) < 1e-4)
            OPTIONS.FreqRange = [];
        end
    end
    
    % ===== INPUT DATA =====
    % If there is nothing special done with the files: files can be handled directly by bst_avg_files
    if isempty(OPTIONS.TimeWindow) && isempty(OPTIONS.ScoutSel) && isempty(OPTIONS.SensorTypes) && isempty(OPTIONS.Rows) && isempty(OPTIONS.FreqRange) && ~OPTIONS.isAvgTime && ~OPTIONS.isAvgRow && ~OPTIONS.isAvgFreq
        InputSetA = {sInputsA.FileName};
        if ~isempty(sInputsB)
            InputSetB = {sInputsB.FileName};
        else
            InputSetB = [];
        end
        OutputType = sInputsA(1).FileType;
    % Else: Call process "Extract values" first
    else
        % Do not concatenate the output
        OPTIONS.Dim = 0;
        % Call extraction process
        [InputSetA, OutputType] = process_extract_values('Extract', sProcess, sInputsA, OPTIONS);
        if ~isempty(sInputsB)
            [InputSetB, tmp__] = process_extract_values('Extract', sProcess, sInputsB, OPTIONS);
        else
            InputSetB = [];
        end
    end

    % === COMPUTE TEST ===
    % Branch between dependend(=paired) and independent tests
    switch (TestType)
        
        % ===== INDEPENDENT TESTS =====
        case {'ttest_equal', 'ttest_unequal', 'abstest', 'abstest2', 'power', 'power_unconstr'}
            % Compute mean and var for both files sets
            [StatA, MessagesA] = bst_avg_files(InputSetA, [], AvgFunction, isAvgVariance, isAvgWeighted, OPTIONS.isMatchRows);
            [StatB, MessagesB] = bst_avg_files(InputSetB, [], AvgFunction, isAvgVariance, isAvgWeighted, OPTIONS.isMatchRows);
            % Add messages to report
            if ~isempty(MessagesA)
                bst_report('Error', sProcess, sInputsA, MessagesA);
                return;
            end
            if ~isempty(MessagesB)
                bst_report('Error', sProcess, sInputsB, MessagesB);
                return;
            end
            if ~isequal(size(StatA.mean), size(StatB.mean))
                bst_report('Error', sProcess, [], 'Files A and B do not have the same number of signals or time samples.');
                return;
            end
            % Do not allow unconstrained sources without a norm
            if (StatA.nComponents ~= 1) && ~OPTIONS.isAbsolute
                bst_report('Error', sProcess, [], ['Cannot run this test on unconstrained sources:' 10 'you must compute the norm of the three orientations first.']);
                return;
            end
            % Bad channels: For recordings, keep only the channels that are good in BOTH A and B sets
            switch lower(sInputsA(1).FileType)
                case 'data'
                    ChannelFlag = StatA.ChannelFlag;
                    ChannelFlag(StatB.ChannelFlag == -1) = -1;
                    isGood = (ChannelFlag == 1);
                case {'results', 'timefreq', 'matrix'}
                    ChannelFlag = [];
                    isGood = true(size(StatA.mean, 1), 1);
                    isGood((StatA.nGoodSamples < 2) | (StatB.nGoodSamples < 2)) = 0;
            end

            % === COMPUTE TEST ===
            % Display progress bar
            bst_progress('start', 'Processes', 'Computing t-test...');
            % Get average results
            mA = StatA.mean(isGood,:,:);
            mB = StatB.mean(isGood,:,:);
            nA = repmat(StatA.nGoodSamples(isGood,:,:), 1, size(mA,2), size(mA,3));
            nB = repmat(StatB.nGoodSamples(isGood,:,:), 1, size(mB,2), size(mB,3));
            % Get variance (if needed)
            if isAvgVariance
                vA = StatA.var(isGood,:,:);
                vB = StatB.var(isGood,:,:);
                % Remove null variances
                iNull = find((vA == 0) | (vB == 0));
                vA(iNull) = eps;
                vB(iNull) = eps;
            else
                iNull = [];
            end
            
            % Compute test statistic
            switch (TestType)
                % === T-TEST: EQUAL VARIANCE ===
                case 'ttest_equal'
                    df   = nA + nB - 2 ;
                    pvar = ((nA-1).*vA + (nB-1).*vB) ./ df;
                    tmap = (mA-mB) ./ sqrt(pvar .* (1./nA + 1./nB));
                    % Calculate p-values from t-values
                    pmap = ComputePvalues(tmap, df, 't', TestTail);
                    % Units: t
                    DisplayUnits = 't';

                % === T-TEST: UNEQUAL VARIANCE ===
                case 'ttest_unequal'
                    df = (vA./nA + vB./nB).^2 ./ ...
                         ((vA./nA).^2./(nA-1) + (vB./nB).^2./(nB-1));
                    tmap = (mA-mB) ./ sqrt(vA./nA + vB./nB);
                    % Calculate p-values from t-values
                    pmap = ComputePvalues(tmap, df, 't', TestTail);
                    % Units: t
                    DisplayUnits = 't';
                    
                % === ABSOLUTE TEST ===
                case 'abstest'
                    % Testing: (|mean(A)| - |mean(B)| = 0)
                    % 
                    % EXPLANATIONS:
                    %   Assume the individual samples xi follow a normal distribution with zero mean and s^2 variance: N(0,s^2).
                    %   Then mean(x) is also normal with zero mean and variance sm^2 = s^2/N.
                    %   
                    %   When we apply abs(mean(x)), we are folding this normal distribution around zero to make it positive. 
                    %   Details are discussed here: https://en.wikipedia.org/wiki/Half-normal_distribution
                    %   If y=|x|, with x~N(0,sm^2), then y has a new distribution with mean my = sm*sqrt(2/pi), and variance sy^2 = sm^2*(1-2/pi)
                    %
                    % RESTRICTIONS
                    %   - A and B are normally distributed (same as t-test assumptions)
                    %   - Cannot be applied if an absolute has been applied already, we need the original values
                    %   
                    % TODO: 
                    % - A and B must be with zero mean, N(0,s^2)  => WRONG
                    %   => Asked Dimitrios, waiting for answer

                    % Test to check that there was no abs already applied, we need the relative values
                    if all(mA(:) > 0) && all(mB(:) > 0)
                        bst_report('Error', sProcess, [], ['This test is designed for values that are positive and negative.' 10 'It cannot be applied to values for which we have already discarded the sign.' 10 'If all your measures you are testing are always strictly positive, then use a Student t-test.']);
                        return;
                    end
            
                    % Mean of: abs(mean(A))-abs(mean(B))
                    mAB = (sqrt(vA./nA) - sqrt(vB./nB)) .* sqrt(2/pi);
                    % Variance of: abs(mean(A))-abs(mean(B))
                    sdAB = sqrt(vA./nA + vB./nB) .* sqrt(1-2/pi) ;
                    S = (abs(mA) - abs(mB) - mAB) ./ sdAB;  %S should be zero mean, unit variance under the null hypothesis (that x and y have zero mean)

                    % [H,P] = ztest(S,0,1);   m = 0; sigma = 1;
                    % zval = (S - m) ./ (sigma ./ sqrt(length(S)));
                    tmap = S .* sqrt(length(S));
                    % Two-tailed test
                    pmap = 2 * (1/2 * erfc(-1 * -abs(tmap) / sqrt(2)));     % 2 * normcdf(-abs(zval),0,1);
                    df = [];
                    % Units: z
                    DisplayUnits = 'z';
                    
                % === ABSOLUTE TEST ===
                case 'abstest2'
                    % Testing: (|mean(A)| - |mean(B)| = 0)
                    % 
                    % EXPLANATIONS:
                    %   Assume the individual samples xi follow a normal distribution with mean "m" and variance "s^2":  X ~ N(m,s^2)
                    %   Then mean(x) is also normal with mean m and variance sm^2 = s^2/N:    mean(X) ~ N(m,s^2/N)
                    %   
                    %   When we apply abs(mean(x)), we are folding this normal distribution to make it positive. 
                    %   Details are discussed here: https://en.wikipedia.org/wiki/Folded_normal_distribution
                    %   If y=|x|, with x~N(m,sm^2), then y has a new distribution with mean my and variance sy^2:
                    %      my = sm*sqrt(2/pi)*exp(-m^2/(2*sm^2)) - m*erf(-m/sqrt(2*sm^2))
                    %         = s/sqrt(N)*sqrt(2/pi)*exp(-m^2/(2*s^2/N)) - m*erf(-m/sqrt(2*s^2/N))
                    %      sy^2 = m^2 + sm^2 - my^2
                    %           = m^2 + s^2/N - my^2
                    %
                    % RESTRICTIONS
                    %   - A and B are normally distributed (same as t-test assumptions)
                    %   - Cannot be applied if an absolute has been applied already, we need the original values
                    %   
                    % TODO: 
                    %   - Not sure this is correct (Asked Dimitrios, waiting for answer)
                    
                    % Test to check that there was no abs already applied, we need the relative values
                    if all(mA(:) > 0) && all(mB(:) > 0)
                        bst_report('Error', sProcess, [], ['This test is designed for values that are positive and negative.' 10 'It cannot be applied to values for which we have already discarded the sign.' 10 'If all your measures you are testing are always strictly positive, then use a Student t-test.']);
                        return;
                    end
            
                    % Mean of: abs(mean(A))-abs(mean(B))
                    % mAabs = sA/sqrt(N)*sqrt(2/pi)*exp(-mA^2/(2*sA^2/nA)) - mA*erf(-mA/sqrt(2*sA^2/nA))
                    %       = sqrt(vA./nA.*(2/pi)) .* exp(-mA.^2/(2.*vA./nA)) - mA*erf(-mA./sqrt(2.*vA./nA))
                    mAabs = sqrt(vA./nA.*(2/pi)) .* exp(-mA.^2./(2.*vA./nA)) - mA.*erf(-mA./sqrt(2.*vA./nA));
                    mBabs = sqrt(vB./nB.*(2/pi)) .* exp(-mB.^2./(2.*vB./nB)) - mB.*erf(-mB./sqrt(2.*vB./nB));
                    mAB = mAabs - mBabs;
                    % Variance of: abs(mean(A))-abs(mean(B))
                    vAabs = mA.^2 + vA./nA - mAabs.^2;
                    vBabs = mB.^2 + vB./nB - mBabs.^2;
                    sdAB = sqrt(vAabs + vBabs);
                    S = (abs(mA) - abs(mB) - mAB) ./ sdAB;  %S should be zero mean, unit variance under the null hypothesis

                    % [H,P] = ztest(S,0,1);   m = 0; sigma = 1;
                    % zval = (S - m) ./ (sigma ./ sqrt(length(S)));
                    tmap = S .* sqrt(length(S));
                    % Two-tailed test
                    pmap = 2 * (1/2 * erfc(-1 * -abs(tmap) / sqrt(2)));     % 2 * normcdf(-abs(zval),0,1);
                    % No need to recompute the values on the fly
                    df = [];
                    % Units: z
                    DisplayUnits = 'z';
                    
                    
                case 'power'
                    % ===== POWER TEST (A/B) =====
                    % Tests for power: test if 2 (z-normalized) sources have same or different power
                    %
                    % If you have xi, n normal random variables with zero mean and unit variance N(0,1), then:
                    % X1 = sum_i(xi^2) is chi-square random variable with n degrees of freedom
                    % https://en.wikipedia.org/wiki/Chi-squared_distribution

                    % If X1 and X2 are chi-square random variables with nA and nB degrees of freedom, then
                    % F = (X1/nA) / (X2/nB) is F-distributed with nA numerator degrees of freedom and nB denominator degrees of freedom.
                    % https://en.wikipedia.org/wiki/F-distribution

                    % The F test can be used to test data sets for difference in power. 
                    % For example, if a source has more power in condition A or condition B. 
                    % It can also be used to test for significant activity in orientation free cases.

                    % But all these rely on having xi normalized N(0,1) sources to begin with. 
                    % Thus, if someone were to produce z-transformed or noise-normalized source maps, then you could create the test:
                    % F = (sum(xi^2) / nA) / (sum(yi^2)/nB)
                    % and then compare it against a parametric F-distribution. See for example:
                    %
                    % (?) ZSCORE ONLY: Shouldn't it be normalized across samples than across time??

                    % F statistic
                    tmap = (mA.^2./nA) ./ (mB.^2./nB);
                    % p value for f-distribution
                    % p = cdf('F',F,nA,nB);  
                    pmap = ones(size(tmap));
                    k = ((tmap > 0) & ~isinf(tmap) & (nA > 0) & (nB > 0));                    
                    pmap(k) = 1 - betainc(nB(k)./(nB(k) + nA(k).*tmap(k)), nB(k)./2, nA(k)./2);
                    % No need to recompute the values on the fly
                    df = [];
                    % Units: F
                    DisplayUnits = 'F';
                    
                case 'power_unconstr'
                    % ===== POWER TEST (TWO SAMPLE / UNCONSTRAINED) =====
                    %  ===== NEED SUM(NORM(X).^2) =====
                    % (?) DOES IT MAKE SENSE ?? 
                    % (?) => Asked Dimitrios, waiting for answer
                    A2 = sum(Ax.^2 + Ay.^2 + Az.^2); %power of sources
                    B2 = sum(Bx.^2 + By.^2 + Bz.^2);
                    F = (A2/(nA*3)) / (B2/(nB*3)); %F statistic
                    p = cdf('F', F, nA*3, nB*3);  % p value for f-distribution 

                    
                %     %% ===== POWER TEST (ONE SAMPLE / UNCONSTRAINED SOURCES) =====
                %     %  ===== NEED NORM(X) =====
                %     %test if an orientation-free source is significant
                %     sx = randn; %x direction, important: z-normalized source
                %     sy = randn; %y direction
                %     sz = randn; %z direction
                %     X = sx.^2 + sy.^2 + sz.^2; %X should be distributed with 3 degrees of freedom
                %     p = cdf('chi',X,3); %p value (if you run this script 1000 times, p should be uniformly distributed under the null)
                %     
                %     (?) Not sure this is really useful for anything.....
                %     (?) => Asked Dimitrios, waiting for answer
                    
                    
                otherwise
                    error('Not supported yet');
            end
            % Remove values with null variances
            if ~isempty(iNull)
                tmap(iNull) = 0;
                pmap(iNull) = 1;
            end
            
            
        % ===== PAIRED/ONE-SAMPLE TESTS =====
        case {'ttest_paired', 'ttest_abspaired', 'ttest_onesample', 'chi2_onesample'}
            % Number of samples must be equal
            if (length(sInputsA) ~= length(sInputsB)) && ~ismember(TestType, {'ttest_onesample', 'chi2_onesample'})
                bst_report('Error', sProcess, [], 'For a paired t-test, the number of files must be the same in the two groups.');
                return;
            end
            % Compute the mean and variance of (samples A - samples B)
            [StatA, MessagesA] = bst_avg_files(InputSetA, InputSetB, AvgFunction, isAvgVariance, isAvgWeighted, OPTIONS.isMatchRows);
            % Add messages to report
            if ~isempty(MessagesA)
                bst_report('Error', sProcess, [], MessagesA);
                return;
            end
            % Display progress bar
            bst_progress('start', 'Processes', 'Computing t-test...');
            % Bad channels and other properties
            switch lower(sInputsA(1).FileType)
                case {'data', 'pdata'}
                    ChannelFlag = StatA.ChannelFlag;
                    isGood = (ChannelFlag == 1);
                case {'results', 'timefreq', 'matrix', 'presults', 'ptimefreq', 'pmatrix'}
                    ChannelFlag = [];
                    isGood = true(size(StatA.mean, 1), 1);
                    isGood(StatA.nGoodSamples < 2) = -1;
            end
            
            % === COMPUTE TEST ===
            % Display progress bar
            bst_progress('start', 'Processes', 'Computing t-test...');
            % Get results
            mean_diff = StatA.mean(isGood,:,:);
            n = repmat(StatA.nGoodSamples(isGood,:,:), 1, size(mean_diff,2), size(mean_diff,3));
            % Get variance (if needed)
            if isAvgVariance
                std_diff = sqrt(StatA.var(isGood,:,:));
                % Remove null variances
                iNull = find(std_diff == 0);
                std_diff(iNull) = eps;
            else
                iNull = [];
            end
            
            % Compute test statistic
            switch (TestType)
                case {'ttest_paired', 'ttest_abspaired', 'ttest_onesample'}
                    % Compute t-test
                    tmap = mean_diff ./ std_diff .* sqrt(n);
                    df = n - 1;
                    % Calculate p-values from t-values
                    pmap = ComputePvalues(tmap, df, 't', TestTail);
                    % Units: t
                    DisplayUnits = 't';
                    
                case {'chi2_onesample'}
                    % https://en.wikipedia.org/wiki/Chi-squared_distribution
                    % =>  If Zi~N(0,1) i=1..n  =>  Q=sum(Zi^2) ~ Chi2(n)
                    % Variable "mean_diff" contains RMS(data)=sqrt(sum(data.^2))   
                    % =>  If data is ~N(0,1), for instance t-stat with (n>30)  =>  mean_diff.^2 ~ Chi2(n)
                    tmap = mean_diff .^ 2;
                    df = n;
                    % Calculate p-values from F-values
                    pmap = ComputePvalues(tmap, df, 'chi2', TestTail);
                    % Units: t
                    DisplayUnits = 'chi2';
                    
                otherwise
                    error('Not supported yet');
            end
            % Remove values with null variances
            if ~isempty(iNull)
                tmap(iNull) = 0;
                pmap(iNull) = 1;
            end
    end
    
    % Return full matrices
    if all(isGood)
        tmap_full = tmap;
        pmap_full = pmap;
        df_full   = df;
    else
        tmap_full = zeros(size(StatA.mean));
        tmap_full(isGood,:,:) = tmap;
        if ~isempty(df)
            df_full = zeros(size(StatA.mean));
            df_full(isGood,:,:) = df;
        end
        if ~isempty(pmap)
            pmap_full = ones(size(StatA.mean));
            pmap_full(isGood,:,:) = pmap;
        end
    end
    
    % === OUTPUT STRUCTURE ===
    % Initialize output structure
    sOutput = db_template('statmat');
    sOutput.pmap         = pmap_full;
    sOutput.tmap         = tmap_full;
    sOutput.df           = df_full;
    sOutput.Correction   = 'no';
    sOutput.Type         = OutputType;
    sOutput.ChannelFlag  = ChannelFlag;
    sOutput.Time         = StatA.Time;
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = DisplayUnits;
    sOutput.nComponents  = StatA.nComponents;
    sOutput.GridAtlas    = StatA.GridAtlas;
    sOutput.Freqs        = StatA.Freqs;
    % Row names
    if isfield(StatA, 'RowNames') && ~isempty(StatA.RowNames)
        if strcmpi(OutputType, 'matrix')
            sOutput.Description = StatA.RowNames;
        elseif strcmpi(OutputType, 'timefreq')
            sOutput.RowNames = StatA.RowNames;
        end
    end
end



%% ===== COMPUTE P-VALUES ====
function p = ComputePvalues(t, df, TestDistrib, TestTail)
    % Default: two-tailed tests
    if (nargin < 4) || isempty(TestTail)
        TestTail = 'two';
    end
    % Default: F-distribution
    if (nargin < 3) || isempty(TestDistrib)
        TestDistrib = 'two';
    end
    % Nothing to test
    if strcmpi(TestTail, 'no')
        p = zeros(size(t));
        return;
    end
    
    % Different distributions
    switch lower(TestDistrib)
        % === T-TEST ===
        case 't'
            % Calculate p-values from t-values 
            switch (TestTail)
                case 'one-'
                    % Inferior one-tailed t-test
                    % p = tcdf(t, df);

                    % Equivalent without the statistics toolbox (FieldTrip formula)            
                    p = 0.5 .* ( 1 + sign(t) .* betainc( t.^2 ./ (df + t.^2), 0.5, 0.5.*df ) );

                case 'two'
                    % Two-tailed t-test
                    % p = 2 * (1 - tcdf(abs(t),df));

                    % Equivalent without the statistics toolbox
                    p = betainc( df ./ (df + t .^ 2), df./2, 0.5);

                    % FieldTrip formula
                    % p2 = 1 - betainc( t.^2 ./ (df + t.^2), 0.5, 0.5.*df );

                case 'one+'
                    % Superior one-tailed t-test
                    % p = 1 - tcdf(t, df);

                    % Equivalent without the statistics toolbox (FieldTrip formula)
                    p = 0.5 .* ( 1 - sign(t) .* betainc( t.^2 ./ (df + t.^2), 0.5, 0.5.*df ) );
            end
            
        % === F-TEST ===
        case 'f'
            error('Not supported yet');
            % Calculate p-values from F-values 
            switch (TestTail)
                case 'one-'
                    % Inferior one-tailed F-test
                    % p = fcdf(t, df(1), df(2));                    
                case 'two'

                case 'one+'
                    % Superior one-tailed F-test
                    % p = 1 - fcdf(t, df(1), df(2));
            end
            
        % === CHI2-TEST ===
        case 'chi2'
            % Calculate p-values from Chi2-values 
            %   chi2cdf(x,n) = gamcdf (x, n/2, 2)
            %    gamcdf(x,a,b) = gammainc(x/b, a)
            %   chi2cdf(x,n) = gammainc(x/2, n/2)
            switch (TestTail)
                case 'one-'
                    % Inferior one-tailed Chi2-test  - NOT USEFUL
                    % p = chi2cdf(t,df);
                    % p = gammainc(t./2, df./2);
                    error('Not relevant.');
                case 'two'
                    % Two-tailed Chi2-test  - NOT USEFUL
                    error('Not relevant.');
                case 'one+'
                    % Superior one-tailed Chi2-test
                    % p = 1 - fcdf(t, df(1), df(2));
                    p = 1 - gammainc(t./2, df./2);
            end
    end
end


    
    
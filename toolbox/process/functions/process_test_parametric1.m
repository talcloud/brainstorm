function varargout = process_test_parametric1( varargin )
% PROCESS_TEST_PARAMETRIC1: Parametric one-sample tests (test vs zero).

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
% Authors: Francois Tadel, Dimitrios Pantazis, 2015-2016
macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Parametric test against zero';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Stat1';
    sProcess.SubGroup    = 'Test';
    sProcess.Index       = 701;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/Statistics';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data',  'results',  'timefreq',  'matrix'};
    sProcess.OutputTypes = {'pdata', 'presults', 'ptimefreq', 'pmatrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
    
    % === GENERIC EXTRACT OPTIONS
    % Label
    sProcess.options.extract_title.Comment    = '<B><U>Select data to test</U></B>:';
    sProcess.options.extract_title.Type       = 'label';
    sProcess.options.extract_title.InputTypes = {'data', 'results', 'timefreq', 'matrix'};
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
    sProcess.options.test_type.Comment = {['<B>One-sample Student''s t-test</B> &nbsp;&nbsp;&nbsp;<FONT color="#777777">X~N(m,s)</FONT><BR>' ...
                                           't = mean(X) ./ std(X) .* sqrt(n) &nbsp;&nbsp;&nbsp;&nbsp; <FONT COLOR="#777777">df=n-1</FONT>'], ...
                                          ['<B>One-sample Chi-square test</B> &nbsp;&nbsp;&nbsp; <FONT color="#777777">Zi~N(0,1), i=1..n</FONT><BR>' ...
                                           'Q = sum(|Zi|^2) &nbsp;&nbsp;&nbsp; <FONT color="#777777">Q~Chi2(n)</FONT>']; ...
                                          'ttest_onesample', 'chi2_onesample'};
    sProcess.options.test_type.Type    = 'radio_label';
    sProcess.options.test_type.Value   = 'ttest_onesample';
    % === TAIL FOR THE TEST STATISTIC
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'two';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = process_test_parametric2('FormatComment', sProcess);
end


%% ===== RUN =====
function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    sOutput = process_test_parametric2('Run', sProcess, sInputs, []);
end






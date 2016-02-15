function F = in_fread_lena(sFile, sfid, iEpoch, SamplesBounds, CorrectionFlag)
% IN_FREAD_LENA:  Read a block of recordings from a LENA file
%
% USAGE:  F = in_fread_lena(sFile, sfid, iEpoch, SamplesBounds, CorrectionFlag) 
%         F = in_fread_lena(sFile, sfid, iEpoch, SamplesBounds) : Read all channels
%         F = in_fread_lena(sFile, sfid, iEpoch)                : Read all channels, all the times
%         F = in_fread_lena(sFile, sfid)                        : Read all channels, all the times, from epoch 1
%
% INPUT:
%     - SamplesBounds : Indices of bounds of the time interval to read
%     - sfid          : Pointer to an open file to read data from

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
% Authors: Sylvain Baillet (2004), Francois Tadel (2008-2011)

%% ===== PARSE INPUTS =====
if (nargin < 5) || isempty(CorrectionFlag)
    CorrectionFlag = 0;
end
if (nargin < 4) || isempty(SamplesBounds)
    iTimes = 1:sFile.header.nTimes;
else
    SamplesBounds = SamplesBounds - sFile.prop.samples(1) + 1;
    iTimes = SamplesBounds(1):SamplesBounds(2);
end
if (nargin < 3) || isempty(iEpoch)
    iEpoch = 1;
end

% Build dimensions to import
dim_to_read = cell(1, length(sFile.header.dimNames));
% Epochs dimension
if ~isempty(sFile.header.iDimEpoch) 
    dim_to_read{sFile.header.iDimEpoch} = iEpoch;
end
% Time dimension
dim_to_read{sFile.header.iDimTime} = iTimes - 1;
% Read file
[sFile.header.lena_tree, F] = fastReadLENA(sFile.header.lena_tree, dim_to_read);
% Remove the first dimension if necessary
F = squeeze(F);






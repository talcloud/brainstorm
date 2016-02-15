function [DataMat, ChannelMat] = in_data_fieldtrip(DataFile)
% IN_DATA_FIELDTRIP: Read recordings from FieldTrip structures (timelocked analysis only).
%
% USAGE:  [DataMat, ChannelMat] = in_data_fieldtrip( DataFile )

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
% Authors: Francois Tadel, 2015

% Get format
[fPath, fBase, fExt] = bst_fileparts(DataFile);
% Initialize returned structure
DataMat = db_template('DataMat');
DataMat.Comment  = fBase;
DataMat.Device   = 'FieldTrip';
DataMat.DataType = 'recordings';
DataMat.nAvg     = 1;

% Load structure
ftMat = load(DataFile);
fields = fieldnames(ftMat);
% If the .time field is not directly accessible, try one level down
if ~isfield(ftMat, 'time') && (length(fields) == 1) && isfield(ftMat.(fields{1}), 'time')
    ftMat = ftMat.(fields{1});  
end
% Check all the required fields
if ~isfield(ftMat, 'time') || ~isfield(ftMat, 'avg') || ~isfield(ftMat, 'label')
    error(['This file is not a valid FieldTrip timelocked structure.' 10 'Missing fields: "time", "avg" or "label".']);
end

% Copy accessible fields
DataMat.Time = ftMat.time;
DataMat.F = ftMat.avg;
if isfield(ftMat, 'var')
    DataMat.Std = sqrt(ftMat.var);
end
if isfield(ftMat, 'dof')
    DataMat.nAvg = max(ftMat.dof(:));
end

% No bad channels information
nChannels = length(ftMat.label);
DataMat.ChannelFlag = ones(nChannels, 1);


% Default channel structure
ChannelMat = db_template('channelmat');
ChannelMat.Comment = 'FieldTrip channels';
ChannelMat.Channel = repmat(db_template('channeldesc'), [1, nChannels]);
% For each channel
for i = 1:nChannels
    if (ftMat.label{i}(1) == 'M')
        ChannelMat.Channel(i).Type = 'MEG';
    else
        ChannelMat.Channel(i).Type = 'EEG';
    end
    ChannelMat.Channel(i).Name    = ftMat.label{i};
    ChannelMat.Channel(i).Loc     = [0; 0; 0];
    ChannelMat.Channel(i).Orient  = [];
    ChannelMat.Channel(i).Weight  = 1;
    ChannelMat.Channel(i).Comment = [];
end





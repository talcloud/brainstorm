function [ftHeadModel, HeadModelMat] = out_fieldtrip_headmodel(HeadModelFile, ChannelFile, iChannels)
% OUT_FIELDTRIP_HEADMODEL: Converts a head model file into a FieldTrip structure (see ft_datatype_headmodel).
% 
% USAGE:  [ftHeadModel, HeadModelMat] = out_fieldtrip_headmodel(HeadModelFile, ChannelFile);
%         [ftHeadModel, HeadModelMat] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat);
%
% INPUTS:
%    - HeadModelFile  : Relative path to a head model file available in the database
%    - HeadModelMat   : Brainstorm head model file structure
%    - ChannelFile    : Relative path to a channel file available in the database
%    - ChannelMat     : Brainstorm channel file structure

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
% Authors: Jeremy T. Moreau, Elizabeth Bock, Francois Tadel, 2015

% ===== LOAD INPUTS =====
% Load head model file
if ischar(HeadModelFile)
    HeadModelMat = in_headmodel_bst(HeadModelFile);
elseif isstruct(HeadModelFile)
    HeadModelMat = HeadModelFile;
else
    error('Failed to load head model.');
end
% Load channel file
if ischar(ChannelFile)
    ChannelMat = in_bst_channel(ChannelFile);
elseif isstruct(ChannelFile)
    ChannelMat = ChannelFile;
else
    error('Failed to load channel file.')
end
% Get sensor type
Modality = ChannelMat.Channel(iChannels).Type;
% Get headmodel type
switch (Modality)
    case 'EEG',   HeadModelMethod = HeadModelMat.EEGMethod;
    case 'ECOG',  HeadModelMethod = HeadModelMat.ECOGMethod;
    case 'SEEG',  HeadModelMethod = HeadModelMat.SEEGMethod;
    case {'MEG','MEG MAG','MEG GRAD'}, HeadModelMethod = HeadModelMat.MEGMethod;
end


% ===== CREATE FIELDTRIP STRUCTURE =====
% Headmodel type
switch (HeadModelMethod)
    case 'meg_sphere'
        ftHeadModel.type = 'singlesphere';
        ftHeadModel.r = HeadModelMat.Param(iChannels(1)).Radii(1);
        ftHeadModel.o = HeadModelMat.Param(iChannels(1)).Center(:)';
    case 'eeg_3sphereberg'
        ftHeadModel.type = 'concentricspheres';
        ftHeadModel.r = HeadModelMat.Param(iChannels(1)).Radii(:)';
        ftHeadModel.o = HeadModelMat.Param(iChannels(1)).Center(:)';
        % Get default conductivities
        BFSProperties = bst_get('BFSProperties');
        ftHeadModel.c = BFSProperties(1:3);
    case 'os_meg'
        ftHeadModel.type = 'localspheres';
        ftHeadModel.r = [HeadModelMat.Param(iChannels).Radii]';
        ftHeadModel.o = [HeadModelMat.Param(iChannels).Center]';
    otherwise
        error('out_fieldtrip_headmodel does not support converting this type of head model.');
end
% Unit and labels
ftHeadModel.unit  = 'm';
ftHeadModel.label = {ChannelMat.Channel(iChannels).Name}';





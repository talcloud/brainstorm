function [CorrVal, CorrSens] = lena_get_ctfcomp(lenaTree, sensor_name)
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
% Authors: Francois Tadel, 2009-2010

% Initialize returned values
CorrVal = [];
CorrSens = [];
% Get node: description / sensor_range / sensor_list
sensor_listId = lena_get_child(lenaTree, root(lenaTree), {'description', 'sensor_range', 'sensor_list'}, 1);

% Get the sensors defintions
sensorsId = children( lenaTree , sensor_listId );
% Find the sensor which name is sensor_name
sSensors = get(lenaTree,children(lenaTree,sensorsId));
iSensor  = find(cellfun(@(c)(strcmp(c.type,'chardata') && strcmp(c.value,sensor_name)), sSensors));
if isempty(iSensor)
    warning(['No description found for sensor "' sensor_name '".']);
    return
end
sensor_uid = sSensors{iSensor}.parent;

% Get the coils correction nodes
sensor_contents = get(lenaTree, children(lenaTree,sensor_uid));
iCorrNodes = find(cellfun(@(c)(strcmp(c.type,'element') && strcmp(c.name,'correction')), sensor_contents));
if isempty(iCorrNodes)
    return
end

% Process all correction nodes
CorrVal  = [];
CorrSens = cell(0);
for i = 1:length(iCorrNodes)
    % Get correction node
    corrNode = sensor_contents{iCorrNodes(i)};
    % Get attribute indice of correction TYPE atttribute
    iAttType = find(cellfun(@(c)strcmp(c.key,  'type'), corrNode.attributes), 1);
    if isempty(iAttType)
        continue;
    end
    % Keep only the 3rd gradient correction values
    if ~ismember(corrNode.attributes{iAttType}.val, {'RB3G', 'G3BR'})
        continue
    end
    % Get elements in correction node
    corrNodeEl = get(lenaTree, corrNode.contents);
    if isempty(corrNodeEl)
        continue
    end  
    % Process all nodes
    for j = 1:length(corrNodeEl)
        % Get the indices of both attributes SENSOR and COEFF
        iAttSensor = find(cellfun(@(c)strcmp(c.key, 'sensor'), corrNodeEl{j}.attributes), 1);
        iAttCoeff  = find(cellfun(@(c)strcmp(c.key, 'coeff'),  corrNodeEl{j}.attributes), 1);
        if isempty(iAttSensor) || isempty(iAttCoeff)
            continue;
        end
        % Fill the returned values
        CorrVal(j)  = str2double(corrNodeEl{j}.attributes{iAttCoeff}.val);
        CorrSens{j} = corrNodeEl{j}.attributes{iAttSensor}.val;
    end
end





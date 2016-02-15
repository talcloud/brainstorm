function childId = lena_get_child( lenaTree, rootId, childPath, isWarning )
% LENA_GET_CHILD: Get a children with a given name, or a path to a child (recursively).

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
% Authors: Francois Tadel, 2009

% Parse inputs
if (nargin < 4)
    isWarning = 1;
end
childId = [];

% Path or name ?
if iscell(childPath)
    childName = childPath{1};
elseif ischar(childPath)
    childName = childPath;
else
    error('Invalid call');
end

% Get the children ids
child = children( lenaTree , rootId );
% Get the children names
names = get( lenaTree , child , 'name' );
% Find the required child name
childId = child(strcmp(names, childName));
% Warning if node not found
if isempty(childId)
    if isWarning
        warning(['Can t find any "' childName '" node.']);
    end
    return
end

% If other path elements to process:
if iscell(childPath) && (length(childPath) > 1)
    childId = lena_get_child( lenaTree, childId, childPath(2:end));
end



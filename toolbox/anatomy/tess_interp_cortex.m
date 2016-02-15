function Wmat = tess_interp_cortex(destSurf, srcSurf, isSameSubject, nbNeighbors)
% TESS_INTERP_CORTEX: 3D nearest-neighbor interpolation between two cortex surfaces,hemisphere by hemisphere
%
% USAGE:  Wmat = tess_interp_cortex(destSurf, srcSurf, isSameSubject, nbNeighbors)
%         Wmat = tess_interp_cortex(destSurf, srcSurf)                     : nbNeighbors = 8, isSameSubject = 0
%
% INPUT:
%    - srcSurf     : Tesselation structure (Faces,Vertices,VertConn,SCS)
%    - destSurf    : Tesselation structure (Faces,Vertices,VertConn,SCS)
%    - nbNeighbors : Number of nearest neighbors to be considered in the interpolation (default is 8)
%    - isSameSubject       : If 1, the two surfaces are smoothed and realigned with an ICP algorithm before interpolation
%
% OUTPUT:
%    - Wmat : Interpolation matrix
%
% NOTE: NCS field in tesselation structures (Normalized coordinate system)
%    - NCS is a structure that locates some fiducial points: AC,PC,IH
%      Each point is defined by a [3x1] matrix in a field named after the point's name.
%    - If the SCS field is present and correctly filled: the surfaces will be aligned based on those information
%      If not, this step would be just skipped

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
% Authors: Francois Tadel, 2010

%% ===== PARSE INPUTS =====
% Check number of arguments
if (nargin < 2) || (nargin > 4)
    error('Usage: Wmat = tess_interp_cortex(destSurf, srcSurf, isSameSubject, nbNeighbors)');
end
% Argument: destSurf, srcSurf
if ~isstruct(destSurf) || ~isstruct(srcSurf) || ~isfield(destSurf, 'VertConn') || ~isfield(srcSurf, 'VertConn')
    error('Source and destination surfaces must be structures with the following fields: Faces, Vertices, VertConn');
end
% Argument: Number of neighbors for interpolation
if (nargin < 3) || isempty(nbNeighbors)
    nbNeighbors = 8;
end
% Argument: isSameSubject
if (nargin < 4) || isempty(isSameSubject)
    isSameSubject = 0;
end
GRAPHICS = 1;


%% ===== SPLIT HEMISPHERES =====
% Display Waitbar
hw = waitbar(0, 'Splitting hemispheres...');
% Split left and right hemispheres
[rHsrc, lHsrc, isConnected(1)]  = tess_hemisplit(srcSurf);
[rHdest,lHdest, isConnected(2)] = tess_hemisplit(destSurf);
% If the hemispheres are connected, process as if there were only one hemisphere (the right one)
skipICP = 0;
if any(isConnected)
    warning('Hemispheres are connected in source or destination surface. Cannot separate hemispheres...');
    % Consider all the vertices as the right hemispheres, and forget about the left
    rHsrc  = union(rHsrc, lHsrc);
    rHdest = union(rHdest, lHdest);
    lHsrc  = [];
    lHdest = [];
    % ICP cancelled if one is connected and not the other
    if ~all(isConnected) && ~isSameSubject
        warning('One surface is connected, not the other. Skipping ICP...');
        skipICP = 1;
    end
end

%% ===== INTER-SUBJECT REGISTRATION =====
if ~isSameSubject
    % ===== ALIGN IN MNI SPACE =====
    srcSurf.Vertices  = cs_convert(srcSurf,  'scs', 'mni', srcSurf.Vertices);
    destSurf.Vertices = cs_convert(destSurf, 'scs', 'mni', destSurf.Vertices);
    if isempty(srcSurf.Vertices) || isempty(destSurf.Vertices)
        if ~isempty(hw) && ishandle(hw)
            close(hw);
        end
        error('Please compute the MNI transformation for both subjects before running this interpolation.');
    end
    
    % ===== SMOOTH SURFACES =====
    % Estimate smoothing parameter
    % nSmoothSrc  = 1000;
    % nSmoothDest = 1000;
    nSmoothSrc  = round(.07 * length(srcSurf.Vertices));
    nSmoothDest = round(.07 * length(destSurf.Vertices));
    % Smooth surfaces
    srcSurf.Vertices  = tess_smooth(srcSurf.Vertices,  0.2, nSmoothSrc,  srcSurf.VertConn, 1);
    destSurf.Vertices = tess_smooth(destSurf.Vertices, 0.2, nSmoothDest, destSurf.VertConn, 1);
    % Show smoothed surfaces
    if GRAPHICS
        showAlignment(destSurf, srcSurf, 'Check surfaces: Before ICP');
    end
    
    % ===== ICP REGISTRATION =====
    if ~skipICP
        % ICP: Options structure
        Options.Verbose=true;
        Options.Registration='Affine';
        % Register right hemisphere
        hw = waitbar(0, hw, 'Registering right hemisphere...');
        srcSurf.Vertices(rHsrc,:) = ICP_finite(destSurf.Vertices(rHdest,:), srcSurf.Vertices(rHsrc,:), Options);
        % Register left hemisphere
        if ~isempty(lHsrc)
            hw = waitbar(0, hw, 'Registering left hemisphere...');
            srcSurf.Vertices(lHsrc,:) = ICP_finite(destSurf.Vertices(lHdest,:), srcSurf.Vertices(lHsrc,:), Options);
        end
        % Show result of ICP alignment
        if GRAPHICS
            showAlignment(destSurf, srcSurf, 'Check surfaces: After ICP');
        end
    end
elseif GRAPHICS
    showAlignment(destSurf, srcSurf, 'Check surfaces');
end


%% ===== SHEPARDS INTERPOLATION =====
% Allocate interpolation matrix
nDest = size(destSurf.Vertices, 1);
nSrc  = size(srcSurf.Vertices, 1);
Wmat = spalloc(nDest, nSrc, nbNeighbors * nDest);
% Close waitbar (bst_shepards has its own waitbar)
if ~isempty(hw) && ishandle(hw)
    close(hw);
end
% Interpolation for right hemisphere
WmatTmp = bst_shepards(destSurf.Vertices(rHdest,:), srcSurf.Vertices(rHsrc,:), nbNeighbors, 0);
Wmat = copySparse(Wmat, WmatTmp, rHdest(:), rHsrc(:));
% Interpolation for left hemisphere
if ~isempty(lHsrc)
    WmatTmp = bst_shepards(destSurf.Vertices(lHdest,:), srcSurf.Vertices(lHsrc,:), nbNeighbors, 0);
    Wmat = copySparse(Wmat, WmatTmp, lHdest(:), lHsrc(:));
end
% Close waitbar
if ~isempty(hw) && ishandle(hw)
    close(hw);
end

end


%% ===== COPY SPARSE MATRIX =====
% Matlab syntax is usually too slow and takes too much memory...
% This function replace the following expression (A and B being to huge but very sparse matrices):
%    A(i,j) = B;
function A = copySparse(A, B, iA, jA)
    [iB, jB] = find(B);
    indB = sub2ind(size(B), iB, jB);
    indA = sub2ind(size(A), iA(iB), jA(jB));
    A(indA) = B(indB);
end

%% ===== SHOW ALIGNMENT =====
function showAlignment(destSurf, srcSurf, wndTitle)
    % Check if progressbar is visible
    isProgressBar = bst_progress('isVisible');
    txt = bst_progress('text');
    % Show surface 1
    [hFig, iDS, iFig, hPatch] = view_surface_matrix(destSurf.Vertices, destSurf.Faces, .4, [1 0 0]);
    set(hPatch, 'EdgeColor', 'r');
    % Show surface 2
    [hFig, iDS, iFig, hPatch] = view_surface_matrix(srcSurf.Vertices, srcSurf.Faces, .4, [], hFig);
    set(hPatch, 'EdgeColor', [.6 .6 .6]);
    set(hFig, 'Name', wndTitle);
    % Restore progress bar
    if isProgressBar
        bst_progress('show');
        bst_progress('text', txt);
    end
end

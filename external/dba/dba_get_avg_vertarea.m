function nbVert = dba_get_avg_vertarea(mesh_ref, mesh_2)
% 
% Compute nbVert for subsampling mesh#2 according to the average vertices area of mesh ref
% 
% Yohan Attal - HM-TC project 2013

disp = 0;

[FaceArea_ref, VertexArea_ref] = tess_area(mesh_ref);
avg_ref = mean(VertexArea_ref);
std_ref = std(VertexArea_ref);
if disp
    figure('color','w'),
    subplot(131), title('mesh ref')
    hold on
    plot(VertexArea_ref), ylabel('Vertex area')
    plot(ones(1,length(VertexArea_ref))*avg_ref,'r')
    plot(ones(1,length(VertexArea_ref))*(avg_ref+std_ref),'g--')
    plot(ones(1,length(VertexArea_ref))*(avg_ref-std_ref),'g--')
end
 
[FaceArea_2, VertexArea_2] = tess_area(mesh_2);
avg_2 = mean(VertexArea_2);
std_2 = std(VertexArea_2);
if disp
    subplot(132), title('mesh 2')
    hold on
    plot(VertexArea_2), ylabel('Vertex area')
    plot(ones(1,length(VertexArea_2))*avg_2,'r')
    plot(ones(1,length(VertexArea_2))*(avg_2+std_2),'g--')
    plot(ones(1,length(VertexArea_2))*(avg_2-std_2),'g--')
end

if avg_2 > avg_ref
    disp('DBA> Error : The order of mesh are not consistent with for subsampling. Try to invert the meshes entry positions.');
    nbVert = -1;
else
    r = 1; 
    cst = 0.025; % subsampling by steps of 2.5%
    while(avg_2 <= avg_ref)
        r = r-cst; 
        [mesh_2_tmp.Faces, mesh_2_tmp.Vertices] = reducepatch(mesh_2.Faces, mesh_2.Vertices, r);
        % [fv_2_tmp.faces, fv_2_tmp.vertices] = reducepatch_subdiv(mesh_2, r);
        [FaceArea_tmp, VertexArea_tmp] = tess_area(mesh_2_tmp);
        avg_2 = mean(VertexArea_tmp);
        std_2 = std(VertexArea_tmp);  
        % check the variance of the mesh & disp warning if necessary
        if (std_2 > 2*std_ref)
            disp('DBA> Warning : The std of the new vertices area is twice time higher than the mesh ref vertices area.');
        end
    end
    nbVert = length(mesh_2_tmp.Vertices);
end

if disp
    subplot(133), title('mesh 2 new')
    hold on
    plot(VertexArea_tmp), ylabel('Vertex area')
    plot(ones(1,length(VertexArea_tmp))*avg_2,'r')
    plot(ones(1,length(VertexArea_tmp))*(avg_2+std_2),'g--')
    plot(ones(1,length(VertexArea_tmp))*(avg_2-std_2),'g--')
end




function [faces, vertices] = reducepatch_subdiv(TessMat, dsFactor)

newNbVertices       = length(TessMat.Vertices);

% Reduce number of vertices
[NewTessMat.Faces, NewTessMat.Vertices] = reducepatch(TessMat.Faces, TessMat.Vertices, dsFactor * 0.94);
% Find the vertices that were kept by reducepatch
[tmp, I, J] = intersect(TessMat.Vertices, NewTessMat.Vertices, 'rows');
% Save initial number for vertices
nVertLowInit = size(NewTessMat.Vertices,1);

% Progress bar
bst_progress('start', 'Resample surface', 'Analyzing surface...');
% Calulate face areas and perimeter
FaceArea  = tess_area(NewTessMat);
% Vertex connectivity, normals, Curvature
VertConn    = tess_vertconn(NewTessMat.Vertices, NewTessMat.Faces);
[VertNormals, FaceNormals] = tess_normals(NewTessMat.Vertices, NewTessMat.Faces, VertConn);
Curvature   = tess_curvature(NewTessMat.Vertices, VertConn, VertNormals);
% Get center of each face
FaceCenter = (NewTessMat.Vertices(NewTessMat.Faces(:,1),:) + NewTessMat.Vertices(NewTessMat.Faces(:,2),:) + NewTessMat.Vertices(NewTessMat.Faces(:,3),:)) ./ 3;
% Get center of mass of the vertices
SurfCenter = mean(NewTessMat.Vertices, 1);
% Get large faces to subdivide (perimeter or area)
iBigFaces = find((sum(FaceNormals .* bst_bsxfun(@minus, FaceCenter, [0 0 SurfCenter(3)]), 2) > 0.04) & ...  % Faces pointing outwards (normal in the same direction as position vector)
    (sum(Curvature(NewTessMat.Faces) > 0, 2) >= 2) & ...       % Curvature has to be > 0
    (FaceArea > mean(FaceArea) + 1*std(FaceArea)));            % Face area threshold
% If there are not enough points to add (white matter): perform search on all the surface
if (length(iBigFaces) < .75 * (newNbVertices - size(NewTessMat.Vertices,1)))
    iBigFaces = find(FaceArea > mean(FaceArea) + 2.5 * std(FaceArea));
end
% Display message
disp(sprintf('BST> Subdividing %d faces from the %d faces generated by reducepatch.', length(iBigFaces), length(NewTessMat.Faces)));

% figure;
% Loop over each face
iRmFaces = [];
bst_progress('start', 'Resample surface', 'Subdividing large faces...', 1, length(iBigFaces));
for i = 1:length(iBigFaces)
    bst_progress('inc', 1);
    % Get the face and, the positions of its vertices, and the center of the face
    f = NewTessMat.Faces(iBigFaces(i),:);
    v = NewTessMat.Vertices(f,:);
    c = mean(v,1);
    
    % === BOUNDING BOX ===
    % Get maximum distance to consider around the face
    dmax = 1.2 * max(sqrt(sum(bst_bsxfun(@minus, v, c) .^ 2, 2)));
    % Select the vertices of the high-res surface that in a small sphere around the center of the face
    iVertBox = find(sum(bst_bsxfun(@minus, TessMat.Vertices, c) .^ 2, 2) < dmax.^2);
    
    %             % Display selection
    %             cla; hold on;
    %             plot3(TessMat.Vertices(iVertBox,1), TessMat.Vertices(iVertBox,2), TessMat.Vertices(iVertBox,3), '.', 'tag', 'ptri');
    %             plot3(c(1), c(2), c(3), '*y');
    %             patch('Vertices', v, 'Faces', [1 2 3], 'FaceColor', 'r');
    %             axis vis3d equal; drawnow; rotate3d on;
    
    % Get the vertices for the target face in the hi-resolution surface
    s1 = find(I(J == f(1)) == iVertBox);
    s2 = find(I(J == f(2)) == iVertBox);
    s3 = find(I(J == f(3)) == iVertBox);
    % Error?
    if isempty(s1) || isempty(s2) || isempty(s3)
        disp(sprintf('BST> Cannot subdivide big face #%d (box too small), skipping...', i));
        continue;
    end
    % Get a subset of the vertex connectivity matrix
    boxVertConn = TessMat.VertConn(iVertBox, iVertBox);
    
    % === FIND PATH TO COMMON NODES ===
    % Expand areas around all the vertices until they all overlap
    iter_max = 20;
    iter = 1;
    sx = [];
    while isempty(sx) && (iter < iter_max)
        s1 = union(s1, find(any(boxVertConn(s1,:),1)));
        s2 = union(s2, find(any(boxVertConn(s2,:),1)));
        s3 = union(s3, find(any(boxVertConn(s3,:),1)));
        sx = intersect(intersect(s1, s2), s3);
        iter = iter + 1;
    end
    
    % Expand areas around all the vertices until they all overlap
    iter_max = 50;
    iter = 1;
    istop = 0;
    d1 = 0;
    d2 = 0;
    d3 = 0;
    while (istop < 2) && (iter <= iter_max)
        % Grow from vertex #1
        i1 = find(any(boxVertConn(s1,:),1));
        s1 = [s1, setdiff(i1,s1)];
        d1 = [d1, i1*0+iter];
        % Grow from vertex #2
        i2 = find(any(boxVertConn(s2,:),1));
        s2 = [s2, setdiff(i2,s2)];
        d2 = [d2, i2*0+iter];
        % Grow from vertex #1
        i3 = find(any(boxVertConn(s3,:),1));
        s3 = [s3, setdiff(i3,s3)];
        d3 = [d3, i1*0+iter];
        % If all the vertices are in the region: stop immediately
        if (length(s1) == length(iVertBox)) && (length(s2) == length(iVertBox)) && (length(s3) == length(iVertBox))
            istop = 10;
            % Do one more iterations after all the vertices are identified
        elseif (istop > 0) || (all(ismember([s2 s3],s1)) && all(ismember([s1 s3],s2)) && all(ismember([s1 s2],s3)))
            istop = istop + 1;
        else
            iter = iter + 1;
        end
    end
    % If an error occured: skip face
    if (iter > iter_max)
        disp(sprintf('BST> Cannot subdivide big face #%d (more than %d nodes distance), skipping...', i, iter_max));
        continue;
    end
    % Take intersection of the three regions
    [sx,ix,jx] = intersect(s1, s2);
    d1 = d1(ix);
    d2 = d2(jx);
    [sx,ix,jx] = intersect(sx, s3);
    d1 = d1(ix);
    d2 = d2(ix);
    d3 = d3(jx);
    dx = d1 + d2 + d3;
    
    %             delete(findobj(0, 'tag', 'ptri')); plot3(TessMat.Vertices(iVertBox(s1),1), TessMat.Vertices(iVertBox(s1),2), TessMat.Vertices(iVertBox(s1),3), '.g', 'tag', 'ptri');
    %             delete(findobj(0, 'tag', 'ptri')); plot3(TessMat.Vertices(iVertBox(s2),1), TessMat.Vertices(iVertBox(s2),2), TessMat.Vertices(iVertBox(s2),3), '.b', 'tag', 'ptri');
    %             delete(findobj(0, 'tag', 'ptri')); plot3(TessMat.Vertices(iVertBox(s3),1), TessMat.Vertices(iVertBox(s3),2), TessMat.Vertices(iVertBox(s3),3), '.y', 'tag', 'ptri');
    
    % === SELECT VERTICES INSIDE THE FACE ===
    % Convert sx back to full indices list
    sx = iVertBox(sx);
    % Keep only the ones that project on the face INSIDE the triangle
    isInside = bst_intriangle(v(1,:), v(2,:), v(3,:), TessMat.Vertices(sx,:));
    sx = sx(isInside);
    dx = dx(isInside);
    if isempty(sx)
        disp(sprintf('BST> Cannot subdivide big face #%d (no candidate in the triangle), skipping...', i));
        continue;
    end
    % plot3(TessMat.Vertices(sx,1), TessMat.Vertices(sx,2), TessMat.Vertices(sx,3), '.y');
    
    % === SELECT VERTICES INSIDE THE FACE ===
    % Keep the closest path to all the nodes
    pathLength = min(dx);
    sx = sx(dx <= pathLength + 2);
    %             plot3(TessMat.Vertices(sx,1), TessMat.Vertices(sx,2), TessMat.Vertices(sx,3), '.y');
    
    % === KEEP THE MOST CENTRAL LOCATION ===
    % Find the closest to the face center
    [d,imin] = min(sqrt(sum(bst_bsxfun(@minus, TessMat.Vertices(sx,:), c) .^ 2, 2)));
    sx = sx(imin);
    % Make sure it is not already in the destination surface
    if ismember(sx, I)
        disp(sprintf('BST> Cannot subdivide big face #%d (vertex already selected), skipping...', i));
        continue;
    end
    %             plot3(TessMat.Vertices(sx,1), TessMat.Vertices(sx,2), TessMat.Vertices(sx,3), 'og');
    
    % === ADD VERTEX ===
    % Add the vertex to the list of vertices
    NewTessMat.Vertices = [NewTessMat.Vertices; c];
    iVertNew = size(NewTessMat.Vertices,1);
    I(end+1) = sx;
    J(end+1) = iVertNew;
    % Add the three new faces to the new surface
    NewTessMat.Faces = [NewTessMat.Faces; ...
        f(1), f(2), iVertNew; ...
        f(1), iVertNew, f(3); ...
        iVertNew, f(2), f(3)];
    iRmFaces(end+1) = iBigFaces(i);
end
% Verify the unicity of the vertex selection
if (length(I) ~= length(unique(I)))
    disp('BST> Error: The same vertex was selected multiple times in the high-resolution brain.');
    disp('BST> Using the basic reducepatch results instead...');
    % Call reducepatch only
    [NewTessFile, iSurface, I, J] = tess_downsize( TessFile, newNbVertices, 'reducepatch');
    return;
end
% Remove the deleted faces
NewTessMat.Faces(iRmFaces,:) = [];

faces       = NewTessMat.Faces;
vertices    = NewTessMat.Vertices;





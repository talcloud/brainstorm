function [Results, OPTIONS] = bst_inverse_linear_2015(HeadModel,OPTIONS)
% BST_INVERSE_LINEAR_2015: Compute an inverse solution (minimum norm, dipole fitting or beamformer)
% USAGE:  [Results,OPTIONS] = bst_inverse_linear_2015(HeadModel, OPTIONS) : Compute mininum operator
%                   OPTIONS = bst_inverse_linear_2015()                   : Return default options
%
% DESCRIPTION:
%     This program computes the whitened and weighted minimum-norm operator,
%     (the wMNE imaging kernel), which is used to compute whitened
%     and weighted minimum-norm estimates (MNE).
%         (e.g., J=wMNEoperator*B; where B is the unwhitened data).
%     It can also compute the whitened and noise-normalized dynamic
%     statistical parametric mapping (dSPM) inverse operator, and/or the
%     whitened standardized low resolution brain electromagnetic tomography
%     (sLORETA) inverse operator, which are used to compute whitened source
%     activity dSPM and sLORETA maps.
%         (e.g., S_dSPM=dSPMoperator*B; where B is the unwhitened data).
%         (e.g., S_sLORETA=sLORETAoperator*B; where B is the unwhitened data).
%
%     The function was originally written with the goal of providing some of the same
%     functionality of the MNE software written by Matti Hamalainen, but no
%     guarantees are made that it performs all computations in the same
%     exact way. It also provides some functionalities not available in the
%     MNE software.
%
%     In March 2015, this function was completely overhauled by John
%     Mosher to create dipole fitting, beamforming, and min norm images all
%     in the same imaging kernel framework. After the imaging kernel is
%     generated, subsequent Process Operations can be applied to further
%     interpret the results, such as finding the optimal dipole in each
%     time slice.
%
% INPUTS:
%    - HeadModel: Array of Brainstorm head model structures
%         |- Gain       : Forward field matrix for all the channels (unconstrained source orientations)
%         |- GridOrient : Dipole orientation matrix
%    - OPTIONS: structure
%         |- NoiseCovMat        : Noise covariance structure
%         |   |- NoiseCov       : Noise covariance matrix
%         |   |- FourthMoment   : Fourth moment (F^2 * F^2'/n)
%         |   |- nSamples       : Number of times samples used to compute those measures
%         |- DataCovMat         : Data covariance structure
%         |   |- NoiseCov       : Data covariance matrix  (F*F'/n)
%         |   |- FourthMoment   : Fourth moment (F^2 * F^2'/n)
%         |   |- nSamples       : Number of times samples used to compute those measures
%         |- ChannelTypes   : Type of each channel (for each row of the Leadfield and the NoiseCov matrix)
%         |- InverseMethod  : {'minnorm', 'gls', 'lcmv'}
%         |- InverseMeasure :    | minnorm: {'amplitude',  'dspm', 'sloreta', 'performance'}
%         |                      |     gls: {'amplitude',   'fit',     'chi', 'performance'}
%         |                      |    lcmv: {'amplitude', 'power',     'nai', 'performance'}
%         |- SourceOrient   : String or a cell array of strings specifying the type of orientation constraints for each HeadModel (default: 'fixed')
%         |- Loose          :    | Value that weights the source variances of the dipole components defining the tangent space of the cortical surfaces (default: 0.2).
%         |- UseDepth       : Flag to do depth weighting (default: 1).
%         |- WeightExp      :    | Order of the depth weighting. {0=no, 1=full normalization, default=0.8}
%         |- WeightLimit    :    | Maximal amount depth weighting (default: 10).
%         |- NoiseMethod    : {'shrink', 'reg', 'diag', 'none'}
%         |- NoiseReg       :    | NoiseMethod='reg' : Amount of regularization
%         |- SnrMethod      : {'rms', 'fixed', 'estimate'}
%         |- SnrRms         :    | SnrMethod='rms'   : RMS source amplitude, in nAm  (Default=1000)
%         |- SnrFixed       :    | SnrMethod='fixed' : Fixed SNR value (Default=9)
%
% OUTPUTS:
%    - Results : Source structure
%       NEW in 2015: if GridOrient is returned empty, the orientation is assumed to be
%       the original grid orientation (usually normal to the surface, or unconstrained in a volume). If
%       it is returned here, then we are returning an optimized orientation
%       that overrides the original orientation.
%    - OPTIONS : Return the modified options


% NOTES (mostly not updated for 2015, see notes below):
%     - More leadfield matrices can be used: the solution will combine all
%       leadfield matrices appropriately.
%
%     - This leadfield structure allows to combine surface and volume
%       source spaces with and without dipole orientation constraints,and
%       with and without area or volumetric current density computations.
%       If using a single sphere headmodel, the silent radial
%       component could be eliminated using the SVD (e.g., use
%       bst_remove_silent.m).
%       NOTE by JCM 2015: Under what considerations would we mix different
%       source spaces? I don't understand combining, for example,
%       cortically constrained surface dipoles with unconstrained volume
%       dipoles. In my 2015 version, I anticipate the user has picked a
%       consistent source modeling space, e.g. unconstrained volume or
%       fixed orientation surface, as dictated by the head model.
%     
%
%     - NoisCov: This should be computed from the pre-stimulus period for
%       averaged ERF data (e.g., using MNE), or from an empty room recording
%       for unaveraged spontaneous or induced data.
%
%     - Orientation constrains for dipoles (.SourceOrient field)
%         - "fixed"     : Dipoles constrained to point normal to cortical surfaces
%         - "free"      : No constraints, dipoles point in x,y,z directions
%         - "loose"     : Source variances of dipole components pointing
%                         tangentially to the cortical surfaces are
%                         multipled by OPTIONS.Loose
%         - "optimal"   : NEW IN 2015: Optimal orientation is returned
%                         in GridOrients for each dipole
%       => OBSOLETE IN 2015: For dealing with multiple source spaces with different types of orientation constraints use for example,
%          OPTION.SourceOrient{1}='fixed';
%          OPTION.SourceOrient{2}='free';
%          This has to correspond with the HeadModel Structure
%
%     - sLORETA: TODO FIX: Output values are multiplied by 1e12 for display in 
%                Brainstorm (time series and cortical maps).

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
% Authors:  Rey Rene Ramirez, Ph.D, 2010-2012
%           Francois Tadel, 2010-2015
%           John Mosher, 2014-2015

% JCM updater tools while massively rebuilding
% COPY THE LATEST VERSION INTO THE BST TOOLBOX
% [status,result] = system('cp -p /Volumes/Local_Shed_4TB/Users/mosher/Dropbox/matlab_dropbox/bst_inverse_linear_2015.m /BST/brainstorm3/toolbox/inverse')
%
% COMPARE THE CURRENT VERSION WITH WHAT IS IN THE TOOLBOX
% visdiff('/Volumes/Local_Shed_4TB/Users/mosher/Dropbox/matlab_dropbox/bst_inverse_linear_2015.m' ,'/BST/brainstorm3/toolbox/inverse/bst_inverse_linear_2015.m')
%
% COMPARE THE CURRENT VERSION WITH THE ORIGINAL STARTING VERSION
% visdiff('/Volumes/Local_Shed_4TB/Users/mosher/Dropbox/matlab_dropbox/bst_inverse_linear_2015_org_Mar26.m' ,'/Volumes/Local_Shed_4TB/Users/mosher/Dropbox/matlab_dropbox/bst_inverse_linear_2015.m');


%% ===== DEFINE DEFAULT OPTIONS =====
Def_OPTIONS.NoiseCovMat    = [];
Def_OPTIONS.DataCovMat     = [];
Def_OPTIONS.ChannelTypes   = {};
Def_OPTIONS.InverseMethod  = 'minnorm';
Def_OPTIONS.InverseMeasure = 'amplitude';
Def_OPTIONS.SourceOrient   = {'free'};
Def_OPTIONS.Loose          = 0.2;
Def_OPTIONS.UseDepth       = 1;
Def_OPTIONS.WeightExp      = 0.5;
Def_OPTIONS.WeightLimit    = 10;
Def_OPTIONS.NoiseMethod    = 'shrink';
Def_OPTIONS.NoiseReg       = .1;
Def_OPTIONS.SnrMethod      = 'fixed';
Def_OPTIONS.SnrRms         = 1000;
Def_OPTIONS.SnrFixed       = 9;
Def_OPTIONS.FunctionName   = [];

% Return the default options
if (nargin == 0)
    Results = Def_OPTIONS;
    return
end

% Sanity check while we build a new version outside of Brainstorm
fprintf('\n\nBST_INVERSE > You are running inverse version %s\n',mfilename('fullpath'))

% Make the default for all the leadfields
% JCM: Why are there multiple headmodels being passed to this function?
numL = size(HeadModel,2);

%TODO: Understand this:
if numL > 1,
    error('BST_INVERSE > Mosher does not understand why we are passing multiple headmodels')
end

Def_OPTIONS.SourceOrient = repmat(Def_OPTIONS.SourceOrient, [1 numL]);

% Copy default options to OPTIONS structure (do not replace defined values)
OPTIONS = struct_copy_fields(OPTIONS, Def_OPTIONS, 0);


%% ===== CHECK FOR INCONSISTENT VALUES =====

% JCM: Legacy checks, should be okay to retain

if (OPTIONS.UseDepth && strcmpi(OPTIONS.InverseMeasure, 'sloreta'))
    disp('BST_INVERSE > Depth weighting is not necessary when using sLORETA normalization, ignoring option UseDepth=1');
    OPTIONS.UseDepth = 0;
end

if (numL ~= length(OPTIONS.SourceOrient))
    error('BST_INVERSE > The number of elements in the HeadModel structure should equal the length of the cell array OPTIONS.SourceOrient.')
end

if ~isempty(OPTIONS.Loose) && (OPTIONS.Loose>=1 || OPTIONS.Loose<=0)
    error('BST_INVERSE > Loose value should be smaller than 1 and bigger than 0, or empty for no loose orientations.')
end

if (OPTIONS.WeightExp > 1) || (OPTIONS.WeightExp < 0)
    error('BST_INVERSE > WeightExp should be a scalar between 0 and 1')
end

if (OPTIONS.NoiseReg > 1) || (OPTIONS.NoiseReg < 0)
    error('BST_INVERSE > NoiseReg should be a scalar between 0 and 1')
end


%% ===== NOISE COVARIANCE REGULARIZATION =====
% Convenient notation
C_noise = OPTIONS.NoiseCovMat.NoiseCov; % all of the channels requested
Var_noise = diag(C_noise); % Diagonal vector of the noise variances per channel
nChannels = length(Var_noise);

% JCM: March 2015, probably don't need this, but will retain
% Detect if the input noise covariance matrix is or should be diagonal
if (norm(C_noise,'fro') - norm(Var_noise,'fro')) < eps(single(norm(Var_noise,'fro'))),
    % no difference between the full matrix and the diagonal matrix
    disp(['BST_INVERSE > Detected diagonal noise covariance, ignoring option NoiseMethod="' OPTIONS.NoiseMethod '".']);
    OPTIONS.NoiseMethod = 'diag';
end

% Commentary April 2015: There is a tendency to apply generic scale values to the
% channels to balance them with regards to units, such as Volts and Tesla.
% But we also have a problem with gradiometer channels vs magnetometers. So
% the "natural" way is to use the channel variances themselves to balance
% out the differences between modalities. But we don't want to do each
% channel, since a dead channel has (near) zero variance. So instead we
% calculate a common variance for each modality, to get us in the ball
% park. So we initially treat each modality as Independent and Identically
% Distributed (IID) and pre-whiten by this to bring the modalities into
% closer alignment.

% Because the units can be different, we need to first balance the
% different types of arrays. How many unique ones do we have:
Unique_ChannelTypes = unique(OPTIONS.ChannelTypes);

% calculate the average variance per channel as a prior
Prior_IID = zeros(length(Unique_ChannelTypes),1);
Prior_Matrix = zeros(nChannels,1); % initialize as a vector

for i = 1:length(Unique_ChannelTypes),
    ndx = find(strcmp(Unique_ChannelTypes(i),OPTIONS.ChannelTypes)); % which ones
    Prior_IID(i) = mean(Var_noise(ndx)); % mean of this type
    % let's get a bit more sophisticated on calculating this IID value
    % We wouldn't want a few extreme values distorting our mean
    % How many channels are there:
    len_ndx = length(ndx);
    % let's toss out the upper and lower values
    ndx_clip = round(len_ndx/10); % 10%
    Variances_This_Modality = sort(Var_noise(ndx));
    % trim to central values
    Variances_This_Modality = Variances_This_Modality((ndx_clip+1):(end-ndx_clip));
    Prior_IID(i) = mean(Variances_This_Modality); % mean of the middle values
    Prior_Matrix(ndx) = sqrt(Prior_IID(i)); % map to this part of the array
end

%TODO: Test for bizarre case of Prior_IID too small

% build whitener to balance out the different types of channels
iPrior_Matrix = spdiags(1./Prior_Matrix,0,nChannels,nChannels);
% make into a matrix itself
Prior_Matrix = spdiags(Prior_Matrix,0,nChannels,nChannels); % recolor

% disable:
iPrior_Matrix = eye(size(iPrior_Matrix));
Prior_Matrix = eye(size(Prior_Matrix));

% Now we can use this iPrior_Matrix to "pre-whiten" the noise covariance
% matrix, and the Prior_Matrix to recolor it.

% Block Whitened noise covariance
Cw_noise = iPrior_Matrix * C_noise * iPrior_Matrix;

% Now the units imbalance between different subarrays is theoretically
% balanced.


%% =============== Now test for the rank of this noise covariance matrix ========
% Before any additional regularizations that may interfere

% FT added 06-Nov-2014
% Make sure the noise covariance matrix is strictly symmetrical
% (the previous operations may cause rounding errors that make the matrix not exactly symmetrical)
%if ~isequal(C_noise, C_noise')
%   C_noise = tril(C_noise) + (tril(C_noise)' - diag(diag(C_noise)));
%end
% JCM April 2015, per FT idea above, we ensure that noise covariance is purely symmetric, but
% simpler method:
Cw_noise = (Cw_noise + Cw_noise')/2;

% So the block whitened noise covariance matrix is Cw_noise.

% Now build the whitener, with no need to worry about cross terms in this
% version.

% Note,impossible to be complex by above symmetry check
[Un,Sn2] = svd(Cw_noise,'econ');
Sn = sqrt(diag(Sn2)); % singular values
tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
Rank_Noise = sum(Sn > tol);

% Rank_Noise = min(Rank_Noise,round(length(Sn)/2))

fprintf('BST_INVERSE > Before regularization, keeping %.0f noise eigenvalues out of %.0f original set\n',Rank_Noise,length(Sn));
Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);

Cw_noise = Un*diag(Sn.^2)*Un'; % possibly deficient matrix now

% With it truncated, see if we need any additional regularizations

%NoiseCovMat       : Noise covariance structure
%         |   |- NoiseCov       : Noise covariance matrix
%         |   |- FourthMoment   : Fourth moment (F^2 * F^2'/n)
%         |   |- nSamples       : Number of times samples used to compute those measures
switch(OPTIONS.NoiseMethod) % {'shrink', 'reg', 'diag', 'none'}
    
    case 'none'
        % TODO We need to clarify that "none" in Regularization means no
        % regularization was applied to the computed Noise Covariance
        % Matrix
        % Do Nothing to Cw_noise
        
    case 'diag'
        Cw_noise = diag(diag(Cw_noise)); % strip to diagonal
        
    case 'reg'
        % The unit of "Regularize Noise Covariance" is as a percentage of
        % the mean variance of the modality. We have already set this to
        % unity by the Prior whitening performed above. So we can just use
        % it directly here, but to ensure it, we multiply the Reg by the
        % maximum of the eigenvalues
        
        % Ridge Regression:
        Cw_noise = Cw_noise + Sn2(1)*OPTIONS.NoiseReg * eye(nChannels);
        
    case 'shrink'
        % Has the user calculated the Fourth Order moments?
        if ~isfield(OPTIONS.NoiseCovMat,'FourthMoment'),
            error(['BST_INVERSE > Please recalculate Noise Covariance, to include Fourth Order Moments']);
        end
        
        % Need to scale the Fourth Moment for the modalities
        
        % FIX: Fix input to give only the channels in OPTIONS.NoiseCovmat
        % FIX: Fix nSamples to be a scalar
        
        FourthMomentw = (iPrior_Matrix.^2) * OPTIONS.NoiseCovMat.FourthMoment * (iPrior_Matrix.^2);
        % use modified version attached to this function
        [Cw_noise,shrinkage]=cov1para_local(Cw_noise,FourthMomentw,OPTIONS.NoiseCovMat.nSamples(1));
        fprintf('\nShrinkage factor is %f\n\n',shrinkage)
        % we now have the "shrunk" whitened noise covariance
        
    otherwise
        error(['BST_INVERSE > Unknown method: NoiseMethod="' OPTIONS.NoiseMethod '"']);
        
end

% JCM April 2015: So now we have calculated the noise covariance matrix across all
% modalities. The prior Brainstorm versions would then DELETE the cross
% covariances between modalities. I suspect this was to handle the problem
% of disparate units, which should have been fixed by the above. So this
% version of bst_inverse_linear will NOT delete the cross terms.


% we ensure that noise covariance is purely symmetric, but
% simpler method:
Cw_noise = (Cw_noise + Cw_noise')/2;

% Rebuild the whitener to include the above regularizations
[Un,Sn2] = svd(Cw_noise,'econ');
Sn = sqrt(diag(Sn2)); % singular values

fprintf('BST_INVERSE > Still Keeping %.0f regularized noise eigenvalues out of %.0f original set\n',Rank_Noise,length(Sn));
Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);

% build full rotating whiteners. We don't expect dramatic reductions in
% rank here, and it's convenient to rotate back to the original space.
% So note that these matrices may not be of full rank.
iWw_noise = Un*diag(1./Sn)*Un'; % inverse whitener
Ww_noise = Un*diag(Sn)*Un'; % forward recolor.

% so Ww * Ww' = Cw, and iWw*iWw' = pinv(Cw)

% so the entire "whitening" process is to first apply the modality
% whitener, then apply the unitless whitener.
% iWw_noise * iPrior_Matrix * d;

% For convenience in the later modeling, let's combine into one whitener

iW_noise = iWw_noise * iPrior_Matrix;
W_noise = Prior_Matrix * Ww_noise;


%% Data Covariance Manipulation, if there is one

% April 2015: We calculate the data covariance in the noise WHITENED
% framework.

if isempty(OPTIONS.DataCovMat),
    iW_data = [];
    W_data = [];
else
    
    % Convenient notation
    C_data = OPTIONS.DataCovMat.NoiseCov; % all of the channels requested
    Var_data = diag(C_data); % Diagonal vector of the noise variances per channel
    
    % quick error check
    if nChannels ~= length(Var_data),
        error('BST_INVERSE > Data covariance is not the same size as noise covariance, something is wrong')
    end
    
    % Whiten the data covariance using the Noise Covariance information
    Cw_data = iW_noise * C_data * iW_noise';
    
    % Now the units imbalance between different subarrays is theoretically
    % balanced.
    
    %         |- DataCovMat         : Data covariance structure
    %         |   |- NoiseCov       : Data covariance matrix  (F*F'/n)
    %         |   |- FourthMoment   : Fourth moment (F^2 * F^2'/n)
    %         |   |- nSamples       : Number of times samples used to compute those measures
    
    
    % TODO: April 2015 Possibly Different Data Covariance Methods
    % (possibly). For April 2015, let's just use what we think is best. The
    % problem with cov1para is that we really need to recalculate the
    % fourth order moments in the whitened framework. So just handle more
    % simply for now.
    
    % We're basically relying on the noise whitener here to have handled
    % deficient cases. That may not be wise. However, the alternative test
    % is to simply put the data covariance into the noise covariance
    % matrix.
    
    % Ensure symmetric
    Cw_data = (Cw_data + Cw_data')/2;
    
    % So the block whitened data covariance matrix is Cw_data.
    
    % Now build the whitener,
    
    [Ud,Sd2] = svd(Cw_data,'econ');
    Sd = sqrt(diag(Sd2)); % singular values
    tol = length(Sd) * eps(single(Sd(1))); % single precision tolerance of standard deviations
    Rank_Data = sum(Sd > tol);
    
    % Rank_Data = min(Rank_Data,round(length(Sd)/4))
    
    
    fprintf('BST_INVERSE > Keeping %.0f data eigenvalues out of %.0f original set\n',Rank_Data,length(Sd));
    Ud = Ud(:,1:Rank_Data);
    Sd = Sd(1:Rank_Data);
    
    Sd = Sd/Sd(end); % set minimum to 1

    
    % build full rotating whiteners. We don't expect dramatic reductions in
    % rank here, and it's convenient to rotate back to the original space.
    % So note that these matrices may not be of full rank.
    iWw_data = Ud*diag(1./Sd)*Ud'; % inverse data whitener
    Ww_data = Ud*diag(Sd)*Ud'; % forward data recolor.
    
    % so Ww * Ww' = Cw, and iWw*iWw' = pinv(Cw)
    
    % so the entire data "whitening" process is to first apply the noise
    % whitener, then apply the data whitener
    % iWw_data * iW_noise * d
    
    % For convenience in the later modeling, let's combine into one whitener
    
    iW_data = iWw_data * iW_noise; % total data whitening
    W_data = W_noise *  Ww_data; % total recoloring
    
    % But recall that these were calculated in a noise whitened framework,
    % so apply the noise whitener first, before using these data whiteners.
    
end




%% ===== Source Model Assumptions ===============
%
% Calculated separately from the data and noise covariance considerations
% number of sources:
NumDipoles = size(HeadModel.GridLoc,1); % number of dipoles (source points)
% (insensitive to the number of components per dipole)

% orientation of each source:
% Each current dipole has an explicit or implied covariance matrix for its
% moment, C_q. Let W_q be the matrix square root, C_q = W_q * W_q'.

% If the HeadModelType is "surface" then there is an orientation available,
% if it is "volume" then an orientation is not pre-defined in
% HeadModel.GridOrient.

% So if the headmodel is "volume", user should not have been able to select
% "fixed" (normal to cortex), nor "loose". Francois has controlled this
% through the GUI interface.
% TODO: We should check this assumption in this code

Wq = cell(1,NumDipoles); % source unit orientations


% Optional Depth Weighting Scalar (used particularly in Min Norms)
% Called for if OPTIONS.UseDepth is true
% Calculate Depth Weighting Scalar
Alpha2 = ones(1,NumDipoles); % default depth weighting
Alpha = ones(1,NumDipoles); % square root

if OPTIONS.UseDepth
    
    % See eq. 6.2.10 of MNE Manual version 2.7 (Mar 2010).
    % Original code had some backflips to check for instabilities.
    % Here we take a simpler direct approach.
    
    % We are assuming unconstrained (three source directions) per source
    % point. We form the unconstrained norm of each point
    ColNorms2 = sum(HeadModel.Gain .* HeadModel.Gain); % square of each column
    SourceNorms2 = sum(reshape(ColNorms2,3,[]),1); % Frobenius norm squared of each source
    
    % Now Calculate the *non-inverted* value
    Alpha2 = SourceNorms2 .^ OPTIONS.WeightExp; % note not negative weight (yet)
    AlphaMax2 = max(Alpha2); % largest squared source norm
    % The limit is meant to keep the smallest from getting to small
    Alpha2 = max(Alpha2, AlphaMax2 ./ (OPTIONS.WeightLimit^2)); % sets lower bound on source strengths
    
    % Now invert it
    Alpha2 = AlphaMax2 ./ Alpha2; % goes from WeightLimit^2 to 1, depending on the inverse strength of the source
    
    Alpha = sqrt(Alpha2);
    % Thus deep weak sources can be amplified up to WeightLimit times greater
    % strength, relative to the stronger shallower sources.
    
    
end

% So, at this point, we have an orientation matrix defined for every source
% point, and we have a possible weighting scalar to apply to it. Each
% source has (as yet unknown) source variance of sigma2. The total source
% covariance we are modeling is
% Cq = Alpha2 * sigma2 * Wq * Wq'; % where sigma2 is unknown and to be
% estimated from the data.

% So we now define Wq to include the optional alpha weighting, as
% calculated above

% Initialize each source orientation

fprintf('BST_INVERSE > Using ''%s'' surface orientations\n',OPTIONS.SourceOrient{1})

for i = 1:NumDipoles,
    
    switch OPTIONS.SourceOrient{1}
        case 'fixed'
           % fprintf('BST_INVERSE > Using constrained surface orientations\n');
            NumDipoleComponents = 1;
            tmp = HeadModel.GridOrient(i,:)'; % 3 x 1
            Wq{i} = tmp/norm(tmp); % ensure unity
            
        case 'loose'
            % fprintf('BST_INVERSE > Using loose surface orientations\n');
            NumDipoleComponents = 3;
            tmp = HeadModel.GridOrient(i,:)'; % preferred direction
            tmp = tmp/norm(tmp); % ensure unity
            tmp_perp = null(tmp'); % orientations perpedicular to preferred
            Wq{i} = [tmp tmp_perp*OPTIONS.Loose]; % weaken the other directions
            
        case 'free'
            % fprintf('BST_INVERSE > Using unconstrained orientations\n');
            NumDipoleComponents = 3;
            Wq{i} = eye(NumDipoleComponents);
            
        case 'optimal'
            % fprintf('BST_INVERSE > Using orientations optimized by data.\n');
            % Note that optimal is initially three, then we optimize in code
            NumDipoleComponents = 3;
            Wq{i} = eye(NumDipoleComponents);
            
        otherwise
            error('BST_INVERSE > Unknown Source Orientation')
    end
    
    % L2norm of Wq in everycase above is 1, (even for loose, since L2 norm
    % of matrix is largest singular value of the matrix).
    
    Wq{i} = Alpha(i)*Wq{i}; % now weight by desired amplifier
    
    % So L2 norm of Wq is equal to the desired amplifier (if any).
end





%% ===== PROCESSING LEAD FIELD MATRICES, WEIGHTS & ORIENTATIONS =====

% put all covariance priors into one big sparese matrix
WQ = blkdiag(sparse(Wq{1}),Wq{2:end});
% (by making first element sparse, we force Matlab to use efficient sparse
% mex function)

% With the above defined, then the whitened lead field matrix is simply

L = iW_noise * (HeadModel.Gain * WQ);

if ~isempty(iW_data),
    Ld = iW_data * (HeadModel.Gain * WQ);
else
    Ld = [];
end

    
% The model at this point is d = {A}{x} + n, where d and n are whitened,
% and x has a covariance prior of unity times unknown lambda (= sigma q 2)

% Every NumDipoleComponents columns of L is a source, which we call A

% (We could optimize this to only do when needed, but it's not too
% expensive in any event, taking less than two seconds for 15000 dipoles)

% Decompose for each source
tic
clear A
A(1:NumDipoles) = deal(struct('Ua',[],'Sa',[],'Va',[],'Ud',[],'Sd',[],'Vd',[]));

for i = 1:NumDipoles,
    % index to next source in the lead field matrix
    ndx = ((1-NumDipoleComponents):0) + i*NumDipoleComponents;
    [Ua,Sa,Va] = svd(L(:,ndx),'econ');

    Sad = diag(Sa); % strip to vector
    tol = length(Sad) * eps(single(Sad(1))); % single precision tolerance
    Rank_Dipole = sum(Sad > tol);
    switch 'notrim'
        case 'trim'
            A(i).Ua = Ua(:,1:Rank_Dipole);
            A(i).Sa = Sa(1:Rank_Dipole,1:Rank_Dipole);
            A(i).Va = Va(:,1:Rank_Dipole);
        case 'notrim'
            Sad((Rank_Dipole+1):end) = 0; % force to perfect zeros
            A(i).Ua = Ua;
            A(i).Sa = diag(Sad);
            A(i).Va = Va;
    end % switch
    
    if ~isempty(Ld), % repeat for data whitened
        [Ud,Sd,Vd] = svd(Ld(:,ndx),'econ');
        
        Sdd = diag(Sd); % strip to vector
        tol = length(Sdd) * eps(single(Sdd(1))); % single precision tolerance
        Rank_Dipole = sum(Sdd > tol);
        switch 'notrim'
            case 'trim'
                A(i).Ud = Ud(:,1:Rank_Dipole);
                A(i).Sd = Sd(1:Rank_Dipole,1:Rank_Dipole);
                A(i).Vd = Vd(:,1:Rank_Dipole);
            case 'notrim'
                Sdd((Rank_Dipole+1):end) = 0; % force to perfect zeros
                A(i).Ud = Ud;
                A(i).Sd = diag(Sdd);
                A(i).Vd = Vd;
        end % switch
    end
end
fprintf('BST_INVERSE > Time to decompose is %.1f seconds\n',toc);
% So L = [A.Ua]*blkdiag(A.Sa)*blkdiag(A.Va)'; % good for any subset too.

% Now do a global decomposition for setting SNR and doing min norms
% Won't need the expensive VL of the SVD, and UL is relatively small.

[UL,SL2] = svd(L*L');
SL2 = diag(SL2);
SL = sqrt(SL2); % the singular values of the lead field matrix
tol = length(SL)*eps(single(SL(1))); % single precision tolerance
Rank_Leadfield = sum(SL > tol);
fprintf('BST_INVERSE > Rank of leadfield matrix is %.0f of %.0f\n',Rank_Leadfield,length(SL));
% but don't trim to rank, we use full matrix below

if ~isempty(Ld),
    
    [ULd,SLd2] = svd(Ld*Ld');
    SLd2 = diag(SLd2);
    SLd = sqrt(SLd2); % the singular values of the lead field matrix
    tol = length(SLd)*eps(single(SLd(1))); % single precision tolerance
    Rank_Leadfield_d = sum(SLd > tol);
    fprintf('BST_INVERSE > Rank of data whitened leadfield matrix is %.0f of %.0f\n',Rank_Leadfield_d,length(SLd));
    % but don't trim to rank, we use full matrix below
end

%% =========== SNR METHODS ================

% Recall our source covariance prior is Cx = Lamda * I, where we have already
% factored out the dipolar covariance, Cq = alpha2 * Lambda * Wq * Wq';
%
% We use the SNR settings to establish a prior on the source variance,
% which is in turn used as essentially a regularizer in the inverse
% methods.
%
% This may seem backwards, i.e. the signal variance with respect
% to the noise variance establishes the SNR, but in inverse processing, we
% may do this backwards, setting an SNR in order to regularize, depending
% on the selection below.

%

switch (OPTIONS.SnrMethod)
    case 'rms'
        % user had directly specified the variance
        Lambda =  OPTIONS.SnrRms^2;
        SNR = Lambda * SL2(1); % the assumed SNR for the entire leadfield
    case 'fixed'
        % user had given a desired SNR, set lambda of the Grammian to achieve it
        SNR = OPTIONS.SnrFixed;
        % use the L2 norm of the whitened Grammian
        Lambda = SNR/(SL(1)^2); % thus SL2(1)*Lambda = desired SNR.
        
    case 'estimate'
        % We use the L2 norm of the whitened data covariance to set the
        % SNR, which in turn sets an implied global source variance
        if isempty(OPTIONS.DataCovMat),
            error('BST_INVERSE > Need to calculate a Data Covariance in order to Estimate SNR');
        end
        % basic idea is SNR = max(eig(OPTIONS.DataCovMat.NoiseCov)) / max(eig(OPTIONS.NoiseCovMat.NoiseCov));
        % which is (signal + noise) / noise. We use the entire lead field
        % matrix to establish the SNR, using both the data whitened and
        % noise whitened grammians.
        
        Lambda = (SL2(1) - SLd2(1)) / (SL2(1)*SLd2(1));
        SNR = Lambda * SLd2(1) / (1 - Lambda * SLd2(1));
        fprintf('BST_INVERSE > Estimated SNR from the whitened data matrix is %.1f\n',SNR);
        
        
    otherwise
        error(['BST_INVERSE > Not supported yet: NoiseMethod="' OPTIONS.SnrMethod '"']);
end
fprintf('BST_INVERSE > Confirm units\n')
if exist('engunits', 'file')
    [LambdaY,tmp,LambdaU] = engunits(sqrt(Lambda));
    fprintf('BST_INVERSE > Assumed RMS of the sources is %.3g %sA-m\n',LambdaY,LambdaU);
else
    fprintf('BST_INVERSE > Assumed RMS of the sources is %g A-m\n', sqrt(Lambda));
end
fprintf('BST_INVERSE > Assumed SNR is %.1f (%.1f dB)\n',SNR,10*log10(SNR));


%% ===== INVERSE METHODS =====

% Generate first the inversions (current dipole time series)

switch lower(OPTIONS.InverseMethod) % {minnorm, gls, lcmv}
    %         |- InverseMethod  : {'minnorm', 'gls', 'lcmv'}
    %         |- InverseMeasure :    | minnorm: {'amplitude',  'dspm', 'sloreta', 'performance'}
    %         |                      |     gls: {'amplitude',   'fit',     'chi', 'performance'}
    %         |                      |    lcmv: {'amplitude', 'power',     'nai', 'performance'}
    
    
    
    case 'minnorm'
        
        % We set a single lambda to achieve desired source variance. In min
        % norm, all dipoles have the same lambda. We already added an
        % optional amplifier weighting into the covariance prior above.
        
        % So the data covariance model in the MNE is now
        % Cd = Lambda * L * L' + I
        % = Lambda *UL * SL * UL' + I = UL * (LAMBDA*SL + I) * UL'
        % so invert Cd and use for kernel
        % xhat = Lambda * L' * inv(Cd)
        % we reapply all of the covariances to put data back in original
        % space in the last step.
        
        % as distinct from GLS, all dipoles have a common data covariance,
        % but each has a unique noise covariance.
        
        switch OPTIONS.InverseMeasure % {'amplitude',  'dspm', 'sloreta', 'performance'}
            case 'amplitude'
                OPTIONS.FunctionName = 'mn';
                % xhat = Lambda * L' * inv(Cd)
                % Kernel = Lambda * L' * inv(Lambda * L * L' + eye(size(L,1))) * iW_noise;
                                
                Kernel = Lambda * L' * (UL * diag(1./(Lambda * SL2 + 1)) * UL')*iW_noise;
                
            case 'performance', % Neyman-Pearson Performance
                OPTIONS.FunctionName = 'mnp';
                % performance is technically a weighted chi-square sum
                % We calculate instead the simpler "GLS" performance for the WF
                
                % Each source has a unique noise covariance
                % p = Ua' * inv(Wv) * d, which complicates the final kernel
                % Let's do it the slow direct way
                
                % The assumed MNE Data Covariance is Lambda * L * L' + eye
                
                CdMNE = Lambda * (L * L') + eye(size(L,1));
                
                Kernel = zeros(size(L')); % preallocate
                
                [UdMNE,SdMNE] = svd(CdMNE);
                iWdMNE = UdMNE * diag(1./sqrt(diag(SdMNE))) * UdMNE'; % data covariance whitener
                
                for i = 1:NumDipoles,
                    % index to next source in the lead field matrix
                    ndx = ((1-NumDipoleComponents):0) + i*NumDipoleComponents;
                    A = L(:,ndx);
                    Cs = Lambda * A * A'; % signal covariance of this source
                    
                    % NOISE COVARIANCE WAY, TOO SLOW
                    % Now Cd = SUM_All Cs + I, by our prewhitening of the noise
                    % matrix. For each source, then Cv = Cd - Cs;
                    switch 'Data' % {'Noise','Data'}
                        case 'Noise'
                                Cv = CdMNE - Cs;
                                [Uv,Sv] = svd(Cv);
                                iWv = Uv*diag(1./sqrt(diag(Sv)))*Uv'; %TODO Stability check, but this is synthetic data
                                [Ua,Sa] = svd(iWv*A,'econ'); % decompose the noise whitened A
                                % SNR = (Sa(1)^2); % TODO save somewhere
                                Kernel(ndx,:) = Ua'*iWv;
                            
                        case 'Data'
                            % instead, use the CAP (data covariance) approach
                            
                            [Ua,Sa] = svd(iWdMNE*A,'econ');
                            
                            Kernel(ndx,:) = diag(1./sqrt(1 - Lambda * diag(Sa).^2))*Ua'*iWdMNE;
                            
                        otherwise
                            error('Wrong performance method switch')
                    end
                    
                    
                    
                end
                
                Kernel = Kernel * iW_noise; % overall whitener
                    
            case 'dspm'
                OPTIONS.FunctionName = 'dspm';
                % ===== dSPM OPERATOR =====
                % =========== NEEDS REWRITING by JCM ======
                
                
                % xhat = Lambda * L' * inv(Cd)
                % Kernel = Lambda * L' * inv(Lambda * L * L' + eye(size(L,1))) * iW_noise;
                
                Kernel = Lambda * L' * (UL * diag(1./(Lambda * SL2 + 1)) * UL');
  
                dspmdiag = sum(Kernel .^2, 2);
                if (NumDipoleComponents == 1)
                    dspmdiag = sqrt(dspmdiag);
                elseif (NumDipoleComponents==3 || NumDipoleComponents==2)
                    dspmdiag = reshape(dspmdiag, NumDipoleComponents,[]);
                    dspmdiag = sqrt(sum(dspmdiag,1)); % Taking trace and sqrt.
                    dspmdiag = repmat(dspmdiag, [NumDipoleComponents, 1]);
                    dspmdiag = dspmdiag(:);
                end
                Kernel = bst_bsxfun(@rdivide, Kernel, dspmdiag);
               
                
                Kernel = Kernel * iW_noise; % overall whitener
                
                
            case 'sloreta'
                OPTIONS.FunctionName = 'sloreta';
                % ===== sLORETA OPERATOR =====
                % =========== NEEDS REWRITING by JCM ======
                
                if false
                    
                    if (NumDipoleComponents(k) == 1)
                        sloretadiag = sqrt(sum(Kernel(start:endd,:) .* L(:,start:endd)', 2));
                        Kernel(start:endd,:) = bst_bsxfun(@rdivide, Kernel(start:endd,:), sloretadiag);
                    elseif (NumDipoleComponents(k)==3 || NumDipoleComponents(k)==2)
                        for spoint = start:NumDipoleComponents(k):endd
                            R = Kernel(spoint:spoint+NumDipoleComponents(k)-1,:) * L(:,spoint:spoint+NumDipoleComponents(k)-1);
                            SIR = sqrtm(pinv(R));
                            Kernel(spoint:spoint+NumDipoleComponents(k)-1,:) = SIR * Kernel(spoint:spoint+NumDipoleComponents(k)-1,:);
                        end
                    end
                    
                    Kernel = Kernel * iW_noise; % overall whitener
                    
                end
                
            otherwise
                error('Unknown Option Inverse Measure: %s',OPTIONS.InverseMeasure)
                
        end
        
        
        
    case 'gls'
        
        % In generalized least-squares, each dipolar source may have a
        % unique lambda, to achieve the desired SNR. So unlike min norm, we
        % may need to set lambda for each and every dipole.
        
        % However, as a prior, this seems counter-intuitive, since it would
        % mean adjusting the variance of deep sources to be much greater
        % (and therefore less reqularized) than shallower sources.
        
        % So in this version, JCM opted to set one global source variance,
        % as in min norm.
        
        % Note that each Wq may already contain a desired amplifier for the
        % gain matrix of that source. That's okay.
        
        % As distinct from minimum norm, every dipole has the exact same
        % noise covariance assumption.
        
        switch OPTIONS.InverseMeasure % {'amplitude',  'fit', 'chi', 'performance'}
            
            
            case 'amplitude'
                OPTIONS.FunctionName = 'gls';
                Kernel = zeros(size(L')); % preallocate
                % recall we already decomposed everything above
                
                for i = 1:NumDipoles,
                    ndx = ((1-NumDipoleComponents):0) + i*NumDipoleComponents;
                    %Kernel(ndx,:) =  A(i).Va*(Lambda*A(i).Sa)*inv(Lambda*A(i).Sa + I)*A(i).Ua';
                    % note this works even for deficient singular values
                    Kernel(ndx,:) =  A(i).Va*(Lambda*A(i).Sa)*diag(1./(diag(Lambda*(A(i).Sa.^2)) + 1))*A(i).Ua';
                end
                
                % we add final whitener and reapply all of the covariances to put data back in original
                % space in the last step.
                
                Kernel = Kernel * iW_noise; % final noise whitening
                                
            case 'performance', % Neyman-Pearson Performance
                OPTIONS.FunctionName = 'glsp';
                Kernel = zeros(size(L')); % preallocate
                
                for i = 1:NumDipoles,
                    ndx = ((1-NumDipoleComponents):0) + i*NumDipoleComponents;
                    %Kernel(ndx,:) =  A(i).Va*(Lambda*A(i).Sa)*inv(Lambda*A(i).Sa + I)*A(i).Ua';
                    Kernel(ndx,:) =  A(i).Ua';
                end
                
                Kernel = Kernel * iW_noise; % final noise whitening

            case 'fit'
                OPTIONS.FunctionName = 'glsfit';
                error('TODO');
                
            case 'chi'
                OPTIONS.FunctionName = 'glschi';
                error('TODO');

            otherwise
                error('Unknown Option Inverse Measure: %s',OPTIONS.InverseMeasure)
        end
        
        
        
    case 'lcmv'
        if isempty(OPTIONS.DataCovMat),
            error('BST_Inverse > Need Data covariance to run lcmv')
        end
        
        switch OPTIONS.InverseMeasure % {'amplitude', 'power',     'nai', 'performance'}
            
            case 'amplitude'
                OPTIONS.FunctionName = 'lcmv';
                % already precomputed
                
                % use the same Lambda as before
                Kernel = zeros(size(L'));
                
                switch 'reg' % {'reg','unreg'}
                    case 'reg'                        
                    
                        Kernel = Lambda*Ld'*iW_data; % the regularized version
                    
                    case 'unreg'
                        % the non-regularized versoin
                        for i = 1:NumDipoles,
                            ndx = ((1-NumDipoleComponents):0) + i*NumDipoleComponents;
                            % already computed
                            Kernel(ndx,:) =  A(i).Vd * pinv(A(i).Sd) * A(i).Ud';
                        end
                        Kernel = Kernel * iW_data;
                        
                        % Final whitener already added in definition of iW_data
                end
                
                
            case 'performance', % Neyman-Pearson Performance
                OPTIONS.FunctionName = 'lcmvp';
                Kernel = zeros(size(L')); % preallocate
                
                for i = 1:NumDipoles,
                    ndx = ((1-NumDipoleComponents):0) + i*NumDipoleComponents;
                    % already computed
                    Kernel(ndx,:) =  diag(1./sqrt(1 - Lambda * diag(A(i).Sd).^2)) * A(i).Ud';
                    % Kernel(ndx,:) =   A(i).Ud'; % simple test case
                end
                
                Kernel = Kernel * iW_data; % final data whitening
                
            case 'nai', % Neural Activity Index of Van Veen
                OPTIONS.FunctionName = 'lcmvnai';
                error('todo');
   
            case 'power'
                OPTIONS.FunctionName = 'lcmvpow';
                error('todo');
                
            otherwise
                error('Unknown Option Inverse Measure: %s',OPTIONS.InverseMeasure)
        end
        
end




%%%%%%% older code  below %%%%%%%%%%%%
if false
    switch 'old'
        case 'minnorm'
            start = 0;
            Kernel = zeros(size(LW,2),size(LW,1)); % note transpose of LW
            lambda2 = 1 ./ SNR; % for regularizing the MN
            for k = 1:numL,
                start = start + 1;
                endd = start + spl(k) - 1;
                Lk = LW(:,start:endd); % this whitened gain matrix
                Rck = Rc(start:endd,start:endd); % covariance priors
                
                % Setup the Min Norm
                % First, generate the population data covariance
                wCD = LW*LW' + (diag(zeros(size(LW,1),1) + lambda2)); % whitened data covariance
                
                % Decompose
                [Ud,Sd] = svd(wCD);
                
                % Data Whitener
                iWd = Ud*diag(1./sqrt(diag(Sd)))*Ud';
                
                % Now process each source in the gain matrix for it's own inversion
                NumSources = size(Lk,2)/NumDipoleComponents(k); % total number of sources
                
                for i = 1:NumSources
                    ndx = ((1-NumDipoleComponents(k)):0)+i*NumDipoleComponents(k); % next index
                    % SVD the data whitened source
                    [Ua,Sa,Va] = svd(iWd*Lk(:,ndx),0); % svd of this source
                    Sa = diag(Sa);
                    Tolerance_Source = size(Lk,1)*eps(single(Sa(1)));
                    Rank_Source = sum(Sa > Tolerance_Source);
                    if Rank_Source < length(Sa),
                        fprintf('%.0f ',i);
                    end
                    % Trim decomposition
                    Ua = Ua(:,1:Rank_Source);
                    Sa = Sa(1:Rank_Source);
                    Va = Va(:,1:Rank_Source);
                    
                    SNR_Weights_Source = Sa.^2 ./ (1 - Sa.^2); % SNR
                    % now write the pseudoinverse results back into this same matrix
                    
                    switch OPTIONS.InverseMeasure
                        case {'amplitude', 'dspm', 'sloreta'}
                            % The true solution
                            Lk(:,ndx) = iWd * Ua * diag(Sa) * Va' * Rck(ndx,ndx);
                        case 'performance'
                            % The model performance
                            %Lk(:,ndx) = iWd * Ua * diag(sqrt(SNR_Weights_Source)); % model performance
                            Lk(:,ndx) = iWd * Ua; % CHEAT model performance
                        otherwise
                            error('Bad Options String %s',OPTIONS.InverseMethod)
                    end
                end
                
                % Now we have a matrix almost ready for use as an imaging kernel
                % Later, below, the whitener will be added
                
                Kernel(start:endd,:) = Lk';
                start = endd;
            end
            
            
            % === DIPOLE FITTING ===
        case 'gls'
            start = 0;
            Kernel = zeros(size(LW,2),size(LW,1)); % note transpose of LW
            lambda2 = 1 ./ SNR; % for regularizing the GLS
            
            for k = 1:numL
                start = start + 1;
                endd = start + spl(k) - 1;
                Lk = LW(:,start:endd); % this whitened gain matrix
                Rck = Rc(start:endd,start:endd); % covariance priors
                % Now process each source in the gain matrix for it's own inversion
                NumSources = size(Lk,2)/NumDipoleComponents(k); % total number of sources
                
                for i = 1:NumSources
                    ndx = ((1-NumDipoleComponents(k)):0)+i*NumDipoleComponents(k); % next index
                    % SVD the source
                    [Ua,Sa,Va] = svd(Lk(:,ndx),0); % svd of this source
                    Sa = diag(Sa);
                    Tolerance_Source = size(Lk,1)*eps(single(Sa(1)));
                    Rank_Source = sum(Sa > Tolerance_Source);
                    % Trim decomposition
                    Ua = Ua(:,1:Rank_Source);
                    Sa = Sa(1:Rank_Source);
                    Va = Va(:,1:Rank_Source);
                    
                    % calculate the weighted subspace
                    Reg_Weights_Source = Sa.^2 ./ (Sa.^2 + lambda2*Sa(1)^2); % regularizer
                    % now write the pseudoinverse results back into this same matrix
                    switch OPTIONS.InverseMeasure
                        case 'amplitude'
                            % The true pseudo-inverse
                            Lk(:,ndx) = Ua * diag(1./Sa) * Va'*Rck(ndx,ndx);
                        case 'performance'
                            % The model performance
                            % Model may be reduced rank
                            Lk(:,ndx) = 0; % zero out
                            Lk(:,ndx(1:size(Ua,2))) = Ua; % model performance
                        case 'chi'
                            error('Not supported yet.');
                        case 'fit'
                            error('Not supported yet.');
                        case 'glsr'
                            % Regularized
                            Lk(:,ndx) = Ua * diag(Reg_Weights_Source./Sa) * Va'*Rck(ndx,ndx);
                        case 'glsr_p'
                            % Regularized model performance
                            % Model may be reduced rank
                            Lk(:,ndx) = 0; % zero out
                            Lk(:,ndx(1:size(Ua,2))) = Ua * diag(sqrt(Reg_Weights_Source)); % model performance
                        otherwise
                            error('Bad Options String %s',OPTIONS.InverseMethod)
                    end
                end
                
                % Now we have a matrix almost ready for use as an imaging kernel
                % Later, below, the whitener will be added
                
                Kernel(start:endd,:) = Lk';
                start = endd;
            end
            
            % === LCMV BEAMFORMER ===
        case 'lcmv'
            error('Not supported yet.');
    end
end




%% ===== NORMALIZATIONS =====

if false % 2015: we don't need this anymore
    start = 0;
    for k = 1:numL
        start = start + 1;
        endd = start + spl(k) - 1;
        
        
        % ===== LOOSE ORIENTATION: RE-ORIENT COMPONENTS =====
        % WARNING: Changing the orientations of the dipoles from MNE to Brainstorm
        % => Make "Unconstrained" and "loose"/"truncated" models directly comparable
        if ~isempty(Q_Cortex{k})
            % Creating a block diagonal matrix
            N = endd - start + 1;
            Nout = N / NumDipoleComponents(k) * 3;
            iRow = reshape(repmat(reshape(1:Nout,3,[]), NumDipoleComponents(k), 1), 1, []);
            iCol = reshape(repmat(1:N,3,1), 1, []);
            Q_Cortex{k} = sparse(iRow, iCol, Q_Cortex{k}(:));
            % Applying orientations
            Kernel(start:endd,:) = Q_Cortex{k} * Kernel(start:endd,:);
        end
        
        start=endd;
    end
end


%% ===== ASSIGN IMAGING KERNEL =====
% Multiply inverse operator and whitening matrix, so no need to whiten data.
% prelmultipley by covariance priors to put back into original domain

% Now the Kernel has been designed for the number of dipole components.
% Ordinarily, to convert back to three d, we would use
% Kernel = WQ * Kernel, which puts all one-d answers back in three-d. But
% that's not optimal for storage. We do need, however, to account for the
% possible Alpha increase in the storage.


% Key assumption is that the source prior is norm 1, as designed in the
% beginning of this program.

if strcmp(OPTIONS.InverseMeasure,'amplitude'),
    % we need to put orientation and weighting back into solution
    
    if NumDipoleComponents == 3, % we are in three-d,
        Kernel = WQ * Kernel; % put the source prior back into the solution
    elseif NumDipoleComponents == 1,
        Kernel = diag(Alpha)*Kernel;  % put possible alpha weighting back in
    end
    
end

% now orient the dipoles, if requested
GridOrient = [];

if strcmp(OPTIONS.SourceOrient{1},'optimal'),
    fprintf('\nBST_INVERSE > Optimizing the orientation');
    % universal method. Could be faster in some methods, but
    % straightforward here
    GridOrient = zeros(NumDipoles,3); % we will return optimal orientation
    TmpKernel = zeros(size(Kernel,1)/3,size(Kernel,2)); % new kernel
    for i = 1:NumDipoles,
        ndx = (-2:0) + i*3;
        [Uorients,Sorients] = svd(Kernel(ndx,:),'econ');
        GridOrient(i,:) = Uorients(:,1)';
        TmpKernel(i,:) = Uorients(:,1)'*Kernel(ndx,:);
    end
    
    % now we have the optimal unit orientations
    Kernel = TmpKernel;
    NumDipoleComponents = 1;
    fprintf('\n'); % done

end


% Return results structure
Results.ImagingKernel = Kernel;
Results.GridOrient = GridOrient; % either empty or optimized
Results.ImageGridAmp  = [];
Results.Whitener      = iW_noise; % full noise covariance whitener
Results.DataWhitener = iW_data;  % full data covariance whitener
Results.nComponents = NumDipoleComponents;


end


%% ==============================================================================
%  ===== HELPER FUNCTIONS =======================================================
%  ==============================================================================


%% ======== Ledoit Shrinkage
% Modified to use the precalculated stats

function [sNoiseCov,shrinkage]=cov1para_local(NoiseCov,FourthMoment,nSamples)
%COV1PARA
% function [sNoiseCov,shrinkage]=cov1para(NoiseCov,FourthMoment,nSamples)
%
% Based on Ledoit's "cov1para" with some modifications
% x is t x n, returns
% sigma n x n
%
% shrinkage is the final computed shrinkage factor, used to weight the
%  i.i.d. prior vs the sample estimate. If shrinkage is specified, it is
%  used on input; else, it's computed.
%
%

% Original code from
% http://www.ledoit.net/cov1para.m
% Original Ledoit comments:
% function sigma=cov1para(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards one-parameter matrix:
%    all variances are the same
%    all covariances are zero
% if shrink is specified, then this const. is used for shrinkage

% Based on
% http://www.ledoit.net/ole1_abstract.htm
% http://www.ledoit.net/ole1a.pdf (PDF of paper)
%
% A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices
% Olivier Ledoit and Michael Wolf
% Journal of Multivariate Analysis, Volume 88, Issue 2, February 2004, pages 365-411
%
% Abstract
% Many economic problems require a covariance matrix estimator that is not
% only invertible, but also well-conditioned (that is, inverting it does
% not amplify estimation error). For large-dimensional covariance matrices,
% the usual estimator - the sample covariance matrix - is typically not
% well-conditioned and may not even be invertible. This paper introduces an
% estimator that is both well-conditioned and more accurate than the sample
% covariance matrix asymptotically. This estimator is distribution-free and
% has a simple explicit formula that is easy to compute and interpret. It
% is the asymptotically optimal convex combination of the sample covariance
% matrix with the identity matrix. Optimality is meant with respect to a
% quadratic loss function, asymptotically as the number of observations and
% the number of variables go to infinity together. Extensive Monte-Carlo
% confirm that the asymptotic results tend to hold well in finite sample.


% Original Code Header, updated to be now from (2014)
% http://www.econ.uzh.ch/faculty/wolf/publications/cov1para.m.zip
%
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards constant correlation matrix
% if shrink is specified, then this constant is used for shrinkage
%
% The notation follows Ledoit and Wolf (2003, 2004)
% This version 04/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is released under the BSD 2-clause license.
%
% Copyright (c) 2014, Olivier Ledoit and Michael Wolf
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Wolf's site now,
% http://www.econ.uzh.ch/faculty/wolf/publications/cov1para.m.zip
% some differences from original Ledoit that confirm Mosher's
% original re-coding.

% % de-mean returns
% [t,n]=size(x);
% meanx=mean(x);
% x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
% Provided
% NoiseCov=(1/t).*(x'*x);

% compute prior
n=size(NoiseCov,1); % number of channels
meanvar=mean(diag(NoiseCov));
prior=meanvar*eye(n); % Note, should be near identity by our pre-whitening

% what we call p
%y=x.^2;
%phiMat=y'*y/t - NoiseCov.^2;
phiMat = FourthMoment - NoiseCov.^2;
phi=sum(sum(phiMat));

% what we call r is not needed for this shrinkage target

% what we call c
gamma=norm(NoiseCov-prior,'fro')^2;

% compute shrinkage constant
kappa=phi/gamma;
% ensure bounded between zero and one
shrinkage=max(0,min(1,kappa/nSamples));

% compute shrinkage estimator
sNoiseCov=shrinkage*prior+(1-shrinkage)*NoiseCov;

end
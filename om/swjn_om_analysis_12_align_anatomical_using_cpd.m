%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: prepare the workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%
if 1
    fprintf('adding basic ecog tools to path \n');
    %addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/MatlabTools'));
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    cpd_make
    %addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    %addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    %addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end

%%
subject_id='AMC029';
root_dir='/mindhive/evlab/u/Shared/ECoG/SUBJECTS/';
template_child_dir='/IMAGING/alignment/';
headshape_child_dir='/IMAGING/MATLAB/';
template_path=strcat([root_dir,subject_id,template_child_dir,'brain_model_raw.mat']);
headshape_path=strcat([root_dir,subject_id,headshape_child_dir,subject_id,'_brain.mat']);
% load template
template_lr = load(template_path);
% template=template_lr.cortex;
template.pos=template_lr.cortex.vert;
template.tri=template_lr.cortex.tri;
%
sub_brain=load(headshape_path);
headshape.pos=sub_brain.cortex.vert;
headshape.tri=sub_brain.cortex.tri;

% cpd setting
opt.method='rigid'; % use nonrigid registration
opt.beta=2;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.2;       % noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=0;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
opt.max_it=5;         % max number of iterations
opt.tol=1e-7;          % tolerance

% do transformation;
[transform,~] = cpd_register(template.pos,headshape.pos, opt);
M = eye(4,4);
M(1:3,1:3) = transform.R;
M(1:3,4)   = transform.t;
transform.M=M;
% save transformation;
save_path=strcat([root_dir,subject_id,template_child_dir,opt.method,'_iter_',num2str(opt.max_it),'_transform.mat']);
save(save_path,'transform');
fprintf('finished tranformation \n');

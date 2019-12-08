%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AMC029
clear all 
close all 
%% 
cfg = [];
cfg.method = 'cortexhull';
cfg.headshape = '/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/segmentation/surf/lh.pial';
cfg.fshome = '/Applications/freesurfer';
if ~exist(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/','AMC029_hull_lh.mat'])
hull_lh = ft_prepare_mesh(cfg);
save(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/','AMC029_hull_lh.mat'], 'hull_lh');
else 
    load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/','AMC029_hull_lh.mat'], 'hull_lh');
end 

%% 
pial = load('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/alignment/brain_model_raw.mat');
pial.pos=pial.cortex.vert;
pial.tri=pial.cortex.tri;
%
sub_brain=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/MATLAB/','AMC029_brain.mat']);
sub_crtx.pos=sub_brain.cortex.vert;
sub_crtx.tri=sub_brain.cortex.tri;
%
alignment_mat=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/alignment/alignment/','rigid_iter_200_transform.mat']);
transform=alignment_mat.transform;
sub_pos_trans = sub_crtx.pos*transform.R'+repmat(transform.t',[size(sub_crtx.pos,1) 1]);
%sub_pos_trans=transform.Y;
sub_crtx.pos=sub_pos_trans;
elec.elecpos=sub_brain.tala.trielectrodes;
elec.chanpos=sub_brain.tala.trielectrodes;
elec.elecpos = elec.elecpos*transform.R'+repmat(transform.t',[size(elec.elecpos,1) 1]);
elec.chanpos = elec.chanpos*transform.R'+repmat(transform.t',[size(elec.chanpos,1) 1]);
labels= [arrayfun(@(x) sprintf('A%d',x),[1:6]','uniformoutput',false);...
        arrayfun(@(x) sprintf('B%d',x),[7:42]','uniformoutput',false);...
        arrayfun(@(x) sprintf('C%d',x),[43:120]','uniformoutput',false);...
        arrayfun(@(x) sprintf('D%d',x),[121:124]','uniformoutput',false);...
        arrayfun(@(x) sprintf('E%d',x),[125:128]','uniformoutput',false)];

elec.label=labels;
figure;
ft_plot_mesh(pial,'facecolor','brain');
ft_plot_mesh(sub_crtx);
shg
ft_plot_sens(elec,'label','label');
view([-55 10]);
material dull;
lighting gouraud;
camlight;
%% 
figure;
ft_plot_mesh(sub_crtx);
ft_plot_mesh(hull_lh,'facecolor','brain','facealpha',.5);
shg
ft_plot_sens(elec,'label','label');
view([-55 10]);
material dull;
lighting gouraud;
camlight;
%% align_electrodes for the subject 
elec_acpc_fr = elec;
grids = {'A*','B*','C*','D*','E*'};
for g = 1:numel(grids)
    cfg = [];
    cfg.channel = grids{g};
    cfg.keepchannel = 'yes';
    cfg.elec = elec_acpc_fr;
    cfg.method = 'headshape';
    cfg.headshape = hull_lh;
    cfg.warp = 'dykstra2012';
    cfg.feedback = 'yes';
    elec_acpc_fr = ft_electroderealign(cfg);
end
    
%% 
figure;
ft_plot_mesh(pial)
ft_plot_sens(elec_acpc_fr);
view([-55 10]);
material dull;
lighting gouraud;
camlight;

%% 
cfg = [];
cfg.channel = {'*'};
cfg.elec = elec_acpc_fr;
cfg.method = 'headshape';
cfg.headshape = '/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/segmentation/surf/lh.pial';
cfg.warp = 'fsaverage';
cfg.fshome = '/Applications/freesurfer';
elec_fsavg_frs = ft_electroderealign(cfg);

save(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/','AMC029_elec_fsavg.mat'], 'elec_fsavg_frs');
%% 
f=figure
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
ft_plot_mesh(fspial_lh);
ft_plot_sens(elec_fsavg_frs);
view([-90 20]);
material dull;
lighting gouraud;
camlight;


 print(f, '-painters','-djpeg', 'example_elec.jpeg');
save([subjID '_elec_fsavg_frs.mat'], 'elec_fsavg_frs');

%% 

% mri=ft_read_mri('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/segmentation/mri/orig.mgz');
% ft_determine_coordsys(mri);
% cfg           = [];
% cfg.method    = 'interactive';
% cfg.coordsys  = 'acpc';
% mri_acpc = ft_volumerealign(cfg, mri);
% %% 
% fsmri_acpc = ft_read_mri('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/segmentation/mri/brain.mgz'); % on Windows, use 'SubjectUCI29_MR_acpc.nii'
% ft_determine_coordsys(fsmri_acpc);
% cfg           = [];
% cfg.method    = 'interactive';
% cfg.coordsys  = 'acpc';
% fmri_acpc = ft_volumerealign(cfg, fsmri_acpc );
% 
% 
% figure
% ft_plot_ortho(fmri_acpc.anatomy, 'transform', fmri_acpc.transform, 'style', 'intersect');
% ft_plot_mesh(sub_crtx,'facecolor','cortex');
% ft_plot_sens(elec, 'label', 'on', 'fontcolor', 'w');
% view([-90 20]);
% 
% material dull;
% lighting gouraud;
% camlight;
% %% 
% cfg            = [];
% cfg.nonlinear  = 'yes';
% cfg.spmversion = 'spm12';
% cfg.spmmethod  = 'new';
% fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);
% %% 
% elec_mni_frv = elec_acpc_fr;
% elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, elec_acpc_fr.elecpos, 'individual2sn');
% elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, elec_acpc_fr.chanpos, 'individual2sn');
% elec_mni_frv.coordsys = 'mni';
% %% 
% figure
% ft_plot_ortho(mri_acpc.anatomy, 'transform', mri_acpc.transform, 'style', 'intersect');
% ft_plot_mesh(sub_crtx);
% ft_plot_sens(elec, 'label', 'on', 'fontcolor', 'w');
% view([-90 20]);
% material dull;
% lighting gouraud;
% camlight;
% 
% %% 
% 
% atlas = ft_read_atlas('/Users/eghbalhosseiniasl1/MyCodes/fieldtrip/template/atlas/brainweb/brainweb_fuzzy.mat');
% figure
% [parcellation] = ft_datatype_parcellation(atlas)
% 
% %ft_plot_ortho(test_acpc.anatomy, 'transform', test_acpc.transform, 'style', 'intersect');
% 
% 
% cfg = [];
% cfg.roi = elec_mni_frv.chanpos(match_str(elec_mni_frv.label, 'LHH1'),:);
% cfg.atlas = atlas;
% cfg.inputcoord = 'mni';
% cfg.output = 'label';
% labels = ft_volumelookup(cfg, atlas);
% [~, indx] = max(labels.count);
% labels.name(indx)
% 
% %% 
% sub_crtx.pos=cortex.vert;
% sub_crtx.tri=cortex.tri;
% ft_plot_mesh(sub_crtx,'facecolor',[0,0,1]);
% view([-55 10]);
% material dull;
% lighting gouraud;
% camlight;
% 
% 
% 
% sub_brain=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/MATLAB/','AMC029_brain.mat']);
% sub_crtx.pos=sub_brain.cortex.vert;
% sub_crtx.tri=sub_brain.cortex.tri;
% 
% ft_plot_mesh(sub_crtx,'color',[1,0,0]);
% 
% shg
% 
% %% 
%   [mri] = ft_read_mri('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/AMC029_elec_r.nii');
%   
%   cfg = [];
% cfg.method = 'interactive';
% cfg.coordsys = 'acpc';
% mri_acpc = ft_volumerealign(cfg, mri);
% 
% %% 
% clear all 
% 
% pial = ft_read_headshape('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/surf/lh.pial');
% test_acpc = ft_read_mri('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/segmentation/mri/orig/001.mgz'); % on Windows, use 'SubjectUCI29_MR_acpc.nii'
% 
% %cfg = [];
% %cfg.method = 'interactive';
% %cfg.coordsys = 'acpc';
% %mri_acpc = ft_volumerealign(cfg, test_acpc);
% 
% 
% figure
% %ft_plot_ortho(test_acpc.anatomy, 'transform', test_acpc.transform, 'style', 'intersect');
% ft_plot_ortho(mri_acpc.anatomy, 'transform', mri_acpc.transform, 'style', 'intersect');
% %ft_plot_ortho(test_acpc.anatomy, 'style', 'intersect');
% 
% ft_plot_mesh(pial,'facecolor','cortex');
% 
% sub_brain=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/MATLAB/','AMC029_brain.mat']);
% sub_crtx.pos=sub_brain.cortex.vert;
% sub_crtx.tri=sub_brain.cortex.tri;
% ft_plot_mesh(sub_crtx,'facecolor','brain');
% %ft_plot_sens(elec, 'label', 'on', 'fontcolor', 'w');
% view([-90 20]);
% material dull;
% lighting gouraud;
% camlight;
% shg

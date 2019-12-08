%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_nature_comm';

%
%% 
electrode_coordinates= [-66.52975346	-33.394	-9.2554;...
-66.39570651	-32.661	-2.4742;...
-67.47125172	-28.807	-6.5063;...
-70.41393952	-10.092	-6.8729;...
-70.37278646	-15.596	-5.9565;...
-69.92651424	-19.083	-7.7892;...
-70.25038437	-15.229	-7.2394;...
-70.35175978	-15.229	-12.921;...
-68.45168156	-24.404	-9.4387];

elec.elecpos=electrode_coordinates;
elec.chanpos=electrode_coordinates;
elec.label=arrayfun(@num2str,[1:size(electrode_coordinates,1)]');
%% align_electrodes for the subject 

%% 
load('/Users/eghbalhosseiniasl1/MyCodes/fieldtrip/template/anatomy/surface_pial_left.mat');
close all 
f=figure;
aspect_ration=9.32./4.13;
y=600;
set(f,'position',[591 455 aspect_ration*y y]);
f.Units = 'normalized';
ax=axes('position',[.02,.1,.3,.3*aspect_ration]);

%fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
%fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(mesh);
ft_plot_sens(elec,'marker','.','facecolor',[0,0,0],'edgecolor',[0,0,0],'elecsize',15,'edgealpha',0);
view([-90 10]);
material dull;
lighting gouraud;
camlight;
ax.Title.String={'electrode position based on MNI from fieldtrip','fieldtrip: MNI template'};
ax.FontSize=14
%
ax=axes('position',[.34,.1,.3,.3*aspect_ration]);
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/cvs_avg35_inMNI152/surf/lh.pial');
%fspial_lh.coordsys = 'fsaverage';
ft_plot_mesh(fspial_lh);
ft_plot_sens(elec,'marker','.','facecolor',[0,0,0],'edgecolor',[0,0,0],'elecsize',15,'edgealpha',0);
view([-90 10]);
material dull;
lighting gouraud;
camlight;
ax.Title.String={'electrode position based on MNI from freesurfer','freesurfer:cvs avg 35 in MNI152'};
ax.FontSize=14
%
freesurfer_brain = ft_read_headshape('/Applications/freesurfer/subjects/cvs_avg35_inMNI152/surf/lh.pial');
%
load('/Users/eghbalhosseiniasl1/MyCodes/fieldtrip/template/anatomy/surface_pial_left.mat');
%
ax=axes('position',[.67,.1,.3,.3*aspect_ration]);
h=ft_plot_mesh(freesurfer_brain,'facecolor',[1,.5,.5],'facealpha',.5);
h.DisplayName='freesurfer:cvs avg 35 in MNI152'
h1=ft_plot_mesh(mesh);
h1.DisplayName='fieldtrip: MNI template';
shg
legend('show')
view([-90 10]);
material dull;
lighting gouraud;
camlight;
ax.Title.String='alignment difference between surfaces';
ax.FontSize=14
if ~exist(strcat(analysis_path))
    mkdir(strcat(analysis_path))
end

%print(f, '-djpeg', strcat(analysis_path,'/average_anatomicals_presentation.jpeg'));
%% 
%% 
froi_loc='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/fROI_nature_comm/fROIs_porj_lh.nii';
fpial_loc='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/fROI_nature_comm/fROIs_porj_lh.nii';
LangLoc_loc='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/fROI_nature_comm/LangParcels_n220_LH.hdr';

pial_path='/Applications/freesurfer/subjects/cvs_avg35_inMNI152/surf/lh.pial';
mri_path='/Applications/freesurfer/subjects/cvs_avg35_inMNI152/mri/brain.mgz';
morph_froi=ft_read_mri(froi_loc);

fspial_lh = ft_read_mri(fpial_loc);
langloc_lh = ft_read_mri(LangLoc_loc);
mri_lh = ft_read_mri(mri_path);
pial=ft_read_headshape(pial_path);

h=ft_plot_mesh(fspial_lh,'facealpha',1);
view([-90 10]);
material dull;
lighting gouraud;
camlight;
shg

%% 

%% 
cfg=[];
cfg.method='ortho';
cfg.funcolormap='jet'

cfg.funparameter   = 'anatomy';
cfg.maskparameter = cfg.funparameter;
%cfg.surffile=pial_path;
ft_sourceplot(cfg,langloc_lh, mri_lh)


%% 
cfg=[];
cfg.view='lomni';
cfg.olayUnits='z';
cfg.showLabels='y';
cfg.pialOverlay='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/fROI_nature_comm/fROIs_porj_lh.mgh';
plotPialSurf('PT001',cfg);


%%
cd('~/Desktop/ielvis/')
cfg=[];
cfg.view='lomni';
cfg.olayUnits='z';
cfg.pialOverlay='/Users/penfield/handMotorLH.mgh';
>>cfgOut=plotPialSurf('PT001',cfg);


clc;clear;close all;
addpath(genpath('./utils'));
% cam_res=8e-6; % assume a watching camera
cam_res=1000;
wavelength=532e-9;
beamWidth=8e-3;
b=Beam(wavelength,beamWidth,'cam_res',cam_res,...
    'profile','gaussian','profile_sigma',0.5);

b.visProfile("Initial Profile",'on_canvas',1);
% dcmos=2*f*wavelength/dxl;

%% free space propogation
b.reset(); % reset to initial parameters
z=200e-3;
b.prop(z);
b.in_prop=0; % if in_prop=1, vis will give the padded result 
b.visProfile();

%% add aperture
b.reset();
b.visProfile("Initial",'on_canvas',1);
b.in_prop=1;
b.prop(50e-3,b.lens(200e-3));
b.visProfile("After prop",'on_canvas',1);

b.Aperture(2e-3);
b.in_prop=0;
b.visProfile("After Aperture",'on_canvas',1);

%% lens focusing (4f front lens)
b.reset(); 
f=500e-3;
b.in_prop=1;
E_out=b.prop(f,b.lens(f));
b.in_prop=0;
b.visProfile("L1 Focus Plane",'on_canvas',1);

% analyze the simulation diffraction limit
I=abs(E_out(round(b.N/2),:)).^2;
% plot(I);
NA=1.5*b.D/(2*f);
d=0.66*b.lambda/(NA);
fprintf("Theoretical resolution limit %.3f um\n",d*1e3);
p=0.84;
d_simu=b.dx*length(I(I>max(I)*p)); % 0.84 in theory
fprintf("Simulated resolution limit (%.2f max I) %.3f um\n",p,d_simu*1e3);

%% 4f system
b.reset();
b.visProfile('Before 4f','on_canvas',1);
f1=400e-3;
f2=200e-3; % 0.5x
b.in_prop=1;
t_seq={b.lens(f1),-1,b.lens(f2)}; % -1 for free space
z_seq={f1,f2,0};
% prop_seq: Ein->t1->z1->tn->zn->Eout
b.prop_seq(t_seq,z_seq);
b.visProfile('After 4f','on_canvas',1);
%% spatial filtering with grating [still have bug]
b.reset();
b.in_prop=1;
f1=300e-3;
f2=150e-3;

pos=[2699,2560];
r=70;
t_spfilter=b.spfilter(r,pos);
b.in_prop=1;
t_l1=b.lens(f1);
t_l2=b.lens(f2).*b.aperture(b.D/2);

t_grating=b.grating(12,8e-6,1,1,0).*b.aperture(b.D/2);
% figure('Color','White'); imshow(angle(t_grating),[]);
b.visProfile("Before SLM Grating",'on_canvas',1)

t1_seq={t_grating,t_l1}; % SLM->f1->L1->f1
z1_seq={f1,f1};
E_out=b.prop_seq(t1_seq,z1_seq); % first prop seq
figure('Color','White'); imshow(abs(E_out).^2,[]);colormap('jet');colorbar;
% b.interact(t_spfilter); % spatial filtering
% t2_seq={-1,t_l2}; % f2->L2->f2
% z2_seq={f2,f2};
% b.prop_seq(t2_seq,z2_seq); % second prop seq
% 
% b.visProfile("After Grating and Propagation",'on_canvas',1);

%% microlens array
b.reset();
f_mla=44.3e-3;
NL=63;
E_out=b.prop(f_mla,b.lens_array(NL,f_mla));
% b.visProfile("After MLA");
I_out=OpticUtil.retrievePad(abs(E_out).^2,b.sz) ;
figure('Color','White');
imshow(I_out,[]);
%% lens and mla
b.reset();
f=150e-3;
f_mla=44.3e-3;
NL=63;
t1_seq={b.lens(f),b.lens_array(NL,f_mla)};
z1_seq={f+f_mla,f_mla};
b.prop_seq(t1_seq,z1_seq);
b.visProfile()

%% SLM Holography [in dev]
img=squeeze(mean(imread('../data/star.png'),3));
img=imresize(img,b.sz);
img_phase=OpticUtil.GS(img,100);
figure;imshow(img_phase,[]);colormap('gray');colorbar;


%% DMD [in dev]
dmd_alpha=24; % 24 degree
dmd_p=13.68e-6;
tmp=rand(64);
figure;imshow(tmp,[])
phase = OpticUtil.expand_img(2*pi*tmp,10);
dmd_img = OpticUtil.lee_hologram(phase,0.1,0);
dmd_img_pad = OpticUtil.pad_img(dmd_img,b.sz);

t_dmd=b.dmd(dmd_img,dmd_p,dmd_alpha,8);
figure;imshow(abs(t_dmd),[]);colorbar;

%% Vortex beam and interference
b.reset();
b.interact(b.vortex(16));
E_total=b.interfere(b.planewave(0.5,0.5));
% b.visProfile('cmap',colormap('jet'))
I = abs(E_total).^2;

% Plot the intensity distribution
figure;
imagesc(b.x, b.y, I);
axis equal tight off;
colormap(jet(256));
%% interference
b1=Beam(wavelength,beamWidth,'cam_res',cam_res,...
    'profile','gaussian','profile_sigma',0.5);
b2=Beam(wavelength,beamWidth,'cam_res',cam_res,...
    'profile','gaussian','profile_sigma',0.5);
b1.interact(b.planewave(1,0));
b1.interact(b.lens(1000e-3));
b1.interfere(b2.E);
b1.visProfile('cmap',colormap('jet'))

%% See dynamically [in dev]
z0=200e-3;
Es=cell(10,1);
for i=1:10
    Es{i}=b.prop(z0);
end
fig=figure('Color','White','MenuBar','none','ToolBar','none','resize','off');
m=moviein(1000);
for i=1:10
    imshow(abs(Es{i}));
    set(gca, 'Position', [0 0 1 1]);
    set(gca,'visible','off')
    m(i+1)=getframe(fig);
end
clc;clear;close all;
addpath(genpath('./utils'));
cam_res=8e-6; % assume a watching camera
wavelength=532e-9;
beamWidth=8e-3;
b=Beam(wavelength,beamWidth,'cam_res',cam_res,...
    'profile','gaussian','profile_sigma',0.8);

b.visProfile("Initial Profile",1);
% dcmos=2*f*wavelength/dxl;

%% free space propogation
z=100e-3;
b.prop(z);
b.visProfile("Prop Freespace",1);

%% lens focusing (4f front lens)

f=200e-3;
E_out=b.prop(f,b.lens(f));
b.visProfile("L1 Focus Plane",1);
% analyze
I=abs(E_out(round(b.N/2),:)).^2;
% plot(I);
NA=1.5*b.D/(2*f);
d=0.66*b.lambda/(NA);
fprintf("Theoretical resolution limit %.3f um\n",d*1e3);
p=0.84;
d_simu=b.dx*length(I(I>max(I)*p)); % 0.84 in theory
fprintf("Simulated resolution limit (%.2f max I) %.3f um\n",p,d_simu*1e3);

%% 4f back lens
f2=100e-3;
b.prop(f2);
b.visProfile("Before L2",1);
b.interact(b.lens(f2));
b.visProfile("After L2",1);

%% spatial filtering 
spf=b.spfilter(10e-6,[0,0]);
% figure;imshow(spf,[]);
f3=100e-3;
f4=100e-3;
b.visProfile("Before L3");
b.prop(f3,b.lens(f3));
b.visProfile("Before Filter");
b.prop(0,spf);
b.visProfile("After Filter");
b.prop(f4,b.lens(f4));
b.visProfile("After L4");

%% microlens array
b.reset();
f_mla=44.3e-3;
NL=63;
E_out=b.prop(f_mla,b.lens_array(NL,f_mla));
b.visProfile("After MLA");

%% grating
b.reset();
amp=0.87;Tx=1;Ty=0;T=12;dx=8e-6;
t_grating=b.grating(T,dx,amp,Tx,Ty);
b.visProfile("Before SLM Grating")

f=200e-3;
b.interact(t_grating);
b.visProfile("After SLM Grating")
b.prop(200e-3);
b.visProfile()
b.prop(f,b.lens(f));
b.visProfile()
%% SLM (holography)
img=squeeze(mean(imread('../data/star.png'),3));
img=imresize(img,b.sz);
img_phase=OpticUtil.GS(img,100);
figure;imshow(img_phase,[]);colormap('gray');colorbar;
%%
b.reset();
b.interact(img_phase);
% b.Aperture();
f=200e-3;
b.prop(f,b.lens(f));
b.visProfile();
%% 4f lens propogation
f1=50e-3;
f2=100e-3;

%% dmd
dmd_alpha=24; % 24 degree
dmd_p=13.68e-6;
tmp=rand(36); %tmp(100)=1;
figure;imshow(tmp,[])
phase = OpticUtil.expand_img(2*pi*tmp,20);
dmd_img = OpticUtil.lee_hologram(phase,0.1);
% dmd_img_pad = OpticUtil.pad_img(dmd_img,b.sz);

t_dmd=b.dmd(dmd_img,dmd_p,dmd_alpha,8);
figure;imshow(angle(t_dmd),[]);colorbar;
%%
f=200e-3;
% close all;
b.reset();
b.prop(f,t_dmd);
b.prop(f,b.lens(f));
b.visProfile("Lens Focus Plane");
spf=b.spfilter(1e-4,[0,0.00045]);
b.prop(0,spf);
b.visProfile("After SPF")
b.prop(f);
b.interact(b.lens(f));
b.visProfile("After 4f");
%%
hdmd=hadamard(4096);
m=120;
p=reshape(hdmd(m,:),[64,64]);
%%
% b.reset();
% f=200e-3;
% phase = OpticUtil.expand_img(2*pi*p,10);
% dmd_img = OpticUtil.lee_hologram(phase,0.1);
% dmd_img_pad = OpticUtil.pad_img(dmd_img,b.sz);
% b.interact(dmd_img_pad);
% b.prop(f,b.lens(f));
% b.visProfile("Lens Focus Plane");
% spf=b.spfilter(1e-3,[-0.002668,0.002668]);
% b.prop(0,spf);
% b.visProfile("After SPF")
% b.prop(f);
% b.interact(b.lens(f));
% b.visProfile("After 4f");
%% test aperture and vis
b.Aperture(1e-3);
b.visProfile("After Aperture");

%% See dynamically
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
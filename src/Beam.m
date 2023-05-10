classdef Beam < handle

properties
    lambda % wavelength
    D % beam width
    A % amplitude
    Phi % phase
    dx % watching camera
    X % camera coordinate x
    Y % camera coordinate y
    pad_factor
    profile_type
    profile_sigma    
    aperture
end
properties (Access = private)
    canvas
    canvas_t
end

properties(Dependent)
    N % resolution in the camera
    sz
    M % padded N for Fourier computation
    E
    x
    y
    k
    dx_ % higher spatial sample frequency corresponding to padding
    Kz % padded propagation wavevector
end


methods
    function obj=Beam(wavelength,beamWidth,options)
        arguments
            wavelength
            beamWidth
            options.cam_res = 8e-6;
            options.profile = "gaussian"
            options.profile_sigma = 0.8
            options.pad_factor = 5
        end
        obj.dx = options.cam_res;
        obj.D = beamWidth;
        obj.lambda = wavelength;
        
        [obj.X,obj.Y]=meshgrid(obj.x,obj.x);
       
%         obj.Pupil=(obj.X.^2+obj.Y.^2)<=(options.aperture/2)^2;
        obj.profile_sigma=options.profile_sigma;
        obj.profile_type=options.profile;
        obj.A = obj.profile(options.profile,options.profile_sigma);
        obj.Phi = zeros(obj.N);
        obj.pad_factor = options.pad_factor;
        obj.aperture=obj.Aperture(obj.D/2);
        obj.canvas=figure('Color','White','Name',"Propagation Figures","Visible","off");
        obj.canvas_t=tiledlayout(obj.canvas,'flow','TileSpacing','none','Padding','none');
    end

    function N=get.N(obj)
        N=round(obj.D/obj.dx);
    end

    function E=get.E(obj)
        E=obj.A.*exp(1i*obj.Phi);

    end

    function k=get.k(obj)
        k=2*pi/obj.lambda;
    end

    function sz=get.sz(obj)
        sz=[obj.N,obj.N];
    end

    function Kz=get.Kz(obj)
        % Essential: use kx that is extended in XYpad 
        % [Use more spatial frequency in calculation]
        % i.e., consider more angles in angular spectrum decomposition
        
        dfx = 1/(obj.dx*obj.M); % real spatial frequency interval
        cx=floor(obj.M/2);
        if mod(obj.N,2)
            fx = (-cx:cx)*dfx; % use more spatial frequency in cal
        else
            fx = (-cx+1:cx)*dfx;
        end
        fy=fx;
        [FX,FY] = meshgrid(fx, fy);
        KX=2*pi*FX; KY=2*pi*FY;
        % Define transfer function
        Kz=sqrt(obj.k^2-KX.^2-KY.^2);
    end

    function x=get.x(obj)
%         x=linspace(-obj.D/2+obj.dx,obj.D/2,obj.N);
        if mod(obj.N,2)
            x=obj.dx*(-floor(obj.N/2):floor(obj.N/2));
        else
            x=obj.dx*(-obj.N/2+1:obj.N/2);
        end
    end

    function y=get.y(obj)
        y=obj.x;
    end

    function M=get.M(obj)
%         M=2^nextpow2(obj.pad_factor*obj.N); % fourier pad
        M=obj.N*obj.pad_factor;
        if mod(obj.N,2)
            if mod(M,2)
                M=M+1;
            end
        end
    end

    function dx_=get.dx_(obj)
        dx_=obj.D/obj.M;
    end

    function set.pad_factor(obj,val)
        obj.pad_factor=val;
    end

    function reset(obj)
        obj.A = obj.profile(obj.profile_type,obj.profile_sigma);
        obj.Phi = zeros(obj.N);
%         obj.aperture=obj.Aperture(obj.D/2);
        obj.canvas=figure('Color','White','Name',"Propagation Figures","Visible","off");
        obj.canvas_t=tiledlayout(obj.canvas,'flow','TileSpacing','none','Padding','none');
    end

    % Characterization
    function A=profile(obj,profile_type,profile_sigma)
        arguments
            obj
            profile_type
            profile_sigma = 0.8
        end
        if profile_type=="gaussian"
            A=OpticUtil.Gaussian(obj.N,0,profile_sigma);
        else
            A=ones(obj.N);
        end
    end
    
    function a=Aperture(obj,r)
        a=sqrt(obj.X.^2+obj.Y.^2)<=r;
        obj.A=obj.A.*a;
        obj.Phi=obj.Phi.*a;

    end

    % Propogation

    function E_out=prop(obj,z,E_mod)
        if nargin<3
            U_pad=OpticUtil.centerPad(obj.E,obj.pad_factor);
        else
            U_pad=OpticUtil.centerPad(obj.E.*E_mod,obj.pad_factor);
        end
        % prop in free space
        
        U_spec=fftshift(fft2(U_pad));
        U_prop=U_spec.*exp(1j*obj.Kz*z);
        E_out = ifft2(ifftshift(U_prop));
        clear U_pad U_spec U_prop;
%         E_out = OpticUtil.retrievePad(E_out,obj.sz);
        
        obj.A = abs(E_out);
        obj.Phi = angle(E_out);
      
    end

    function E_out=prop_seq(obj,t_seq,z_seq,options)
        % z_seq: prop distance sequence
        % t_seq: complex amplitude sequence
        % seq: Ein->t1->z1->tn->zn->Eout
        arguments
            obj
            t_seq
            z_seq
            options.inputpad = 1
        end
        if options.inputpad
            E_in=OpticUtil.centerPad(obj.E,obj.pad_factor);
        else
            E_in=obj.E;
        end

        for i=1:length(t_seq)
            t=t_seq{i};
            z=z_seq{i};
            if z==-1 % spatial filtering
                U_spec=fftshift(fft2(E_in));
                U_prop=U_spec.*t;
            elseif t==-1 % free space prop
                U_spec=fftshift(fft2(E_in));
                U_prop=U_spec.*exp(1j*obj.Kz*z);
            else % interact and prop
                t_pad=OpticUtil.centerPad(t,obj.pad_factor);
                U_spec=fftshift(fft2(E_in.*t_pad));
                U_prop=U_spec.*exp(1j*obj.Kz*z);
            end
            E_out = ifft2(ifftshift(U_prop));
            E_in = E_out;
        end
%         E_out = OpticUtil.retrievePad(E_out,obj.sz);
        clear U_spec U_prop;
        obj.A = abs(E_out);
        obj.Phi = angle(E_out);
     
    end
    % Beam Interaction
    function E_out=interact(obj,t)
        arguments
           obj
           t
        end
        E_out=obj.E.*t;
        obj.A = abs(E_out);
        obj.Phi = mod(angle(E_out),2*pi);
%         if options.objtype=="spfilter"
%             [h,w]=size(obj.A);
%             hc=options.pos(2);
%             wc=options.pos(1);
%             dn=options.pos(3);
%             
%             tmp=zeros(h,w);
%             tmp(h/2-dn:h/2+dn,w/2-dn:w/2+dn)=...
%                 obj.E(hc-dn:hc+dn,wc-dn:wc+dn);
%             obj.A=abs(tmp);
%             obj.Phi=angle(tmp);
%             
%         end
    end

    function E_out=interfere(obj,t)
        E_out=obj.E+t;
        obj.A = abs(E_out);
        obj.Phi = angle(E_out);
    end
    
    
    
    
    % Common amplitude profile
    function t=spfilter(obj,r,pos,N)
    % spatial filter
    % pos: position in real space [xf,yf]
%         t=(obj.X-pos(1)).^2+(obj.Y-pos(2)).^2<=r^2;
        
        if nargin<4
            x1=1:obj.N*obj.pad_factor;
            
        else
            x1=1:N;
        end
        [X1,Y1]=meshgrid(x1,x1);
        t=((X1-pos(1)).^2+(Y1-pos(2)).^2)<=r^2;
    end

    function t=planewave(obj,ux,uy)
        t=exp(1j*obj.k*(ux*obj.X+uy*obj.Y));
    end

    % Common phase profile
    function t=lens(obj,f)
        t=exp(-1j*obj.k*(obj.X.^2+obj.Y.^2)/(2*f));
    end
    
    function t=lens_array(obj,NL,f)
        % NL: microlens pitch (in pixel)
        microlens_phase = zeros(obj.N);
        if mod(NL,2)==1
            xl=(-floor(NL/2):floor(NL/2))*obj.dx;
        else
            xl = (-floor(NL/2)+1:floor(NL/2))*obj.dx;
        end
        [XL,YL] = meshgrid(xl, xl);
        n=ceil(obj.N/NL);
        dn=mod(obj.N,NL);
        for i = 1:n
            for j = 1:n
                yrange=NL*(i-1)+1:NL*i;
                xrange=NL*(j-1)+1:NL*j;
                if(i==n||j==n)
                    ylrange=1:NL; xlrange=1:NL;
                    if (i==n), yrange=NL*(i-1)+1:NL*(i-1)+dn; ylrange=1:dn; end
                    if (j==n), xrange=NL*(j-1)+1:NL*(j-1)+dn; xlrange=1:dn; end
                    microlens_phase(yrange,xrange) = -obj.k*(XL(ylrange,xlrange).^2+YL(ylrange,xlrange).^2)/(2*f);                
                else
                    microlens_phase(yrange,xrange) = -obj.k*(XL.^2+YL.^2)/(2*f);
                end
            end
        end
%         figure;imshow(microlens_phase,[]);colorbar;
        t=exp(1j*microlens_phase);
%         t=t.*obj.aperture;
    end

    function t=vortex(obj, m)
        
        % Calculate azimuthal angle
        theta = atan2(obj.Y, obj.X);
        t = exp(1i*m*theta);

    end

    function t=grating(obj, T, dx, A, Tx, Ty)
        arguments
            obj
            T
            dx
            A = 1
            Tx = 1
            Ty = 0
            
        end
        % Tx: whether use x direction grating
        % Ty: whether use y direction grating
        
        T = T*dx; % 闪耀光栅周期
        grating_phase_x = 2*pi*mod(Tx*obj.X,T)/T;
        grating_phase_y = 2*pi*mod(Ty*obj.Y,T)/T;
        grating_phase = A*mod(grating_phase_x+grating_phase_y,2*pi);
        t=exp(1j*grating_phase);
%         t=t.*obj.aperture;
    end

    function t=dmd(obj,dmd_img,p,alpha,NL_)
        % alpha: reflection angle (in degree)
%         % TODO: consider the 45degree 
    
        dmd_img_pad=OpticUtil.pad_img(dmd_img,obj.sz);
        sin_alp=sin(alpha/360*2*pi);
        dmd_phase = zeros(obj.N);
        NL=round(NL_*p/obj.dx);
        n=ceil(obj.N/NL);
        dn=mod(obj.N,NL);
        for i = 1:n % y
            for j = 1:n % x
                yrange=NL*(i-1)+1:NL*i;
                xrange=NL*(j-1)+1:NL*j;
                if(i==n||j==n)
                    if (i==n), yrange=NL*(i-1)+1:NL*(i-1)+dn; end
                    if (j==n), xrange=NL*(j-1)+1:NL*(j-1)+dn; end
                end
                dmd_phase(yrange,xrange) = obj.k*p*sin_alp*(i-round(n/2)); %sqrt((i-round(n/2))^2+(j-round(n/2))^2);    
            end
        end
        t=exp(1j*dmd_phase).*dmd_img_pad;
%         t=t.*obj.aperture;
    end

    % Visualization
    function visProfile(obj,figname,options)
        arguments
            obj
            figname = 'Beam Profile';
            options.on_canvas = 0;
            options.cmap = 'gray';
        end
        if options.on_canvas
            nexttile(obj.canvas_t,[1,2]);
            obj.canvas.Visible=1;
            montage({obj.A/max(obj.A,[],'all'),mod(obj.Phi,2*pi)/(2*pi)}, 'Size', [1 2],'DisplayRange', []);
            colormap(options.cmap);colorbar;
            title(strcat(figname,"$$A|\Phi/2\pi$$"),'Interpreter','latex');
           
        else 
            figure('Color','White','Name',figname);
            subplot(121);
            imagesc(obj.x, obj.y, obj.A);title('Amplitude');
            axis equal tight;
            colormap(options.cmap);colorbar;
            
            
            axis square
            subplot(122);
            imagesc(obj.x, obj.y, obj.Phi);title('Phase');
            axis equal tight;
            colormap(options.cmap);colorbar;
        end
        
    end
end


end
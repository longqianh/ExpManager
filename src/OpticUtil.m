classdef OpticUtil
    methods(Static)
        function I=Gaussian(D,mu,sigma)
            % D: beamWidth (in pixel)
            y=linspace(-1,1,D);
            x=linspace(-1,1,D);
            [xx,yy]=meshgrid(x,y);
            I = 1/(2*pi*sigma)*exp(-((xx-mu).^2+(yy-mu).^2)/(2*sigma^2));
            I_max=max(max(I));
            I_min=min(min(I));
            I = (xx.^2+yy.^2<=1).*(I-I_min)/(I_max-I_min);
        end

        function I_pad=fourierPad(I,pad_factor)
            src_sz=size(I);
            pad_sz=2.^nextpow2(src_sz*pad_factor);
            I_pad=zeros(pad_sz);
            I_pad(floor(pad_sz(1)/2-src_sz(1)/2)+1:floor(pad_sz(1)/2+src_sz(1)/2),...
                floor(pad_sz(2)/2-src_sz(2)/2)+1:floor(pad_sz(2)/2+src_sz(2)/2))=I;
        end

        function I=retreivePad(I_pad,src_sz)
            pad_sz=size(I_pad);
            I=I_pad(floor(pad_sz(1)/2-src_sz(1)/2)+1:floor(pad_sz(1)/2+src_sz(1)/2),...
                floor(pad_sz(2)/2-src_sz(2)/2)+1:floor(pad_sz(2)/2+src_sz(2)/2));
        end


        function [Uout,xout,yout] = prop_Angular_Spectrum(Uin,lambda,z,dx)
            % 角谱传输仿真函数
            % 输入：
            % Uin: 输入场（2D矩阵）
            % lambda: 波长
            % z: 传播距离
            % dx: 输入场的空间采样间隔
            % 输出：
            % Uout: 传播后的场（2D矩阵）
            % xout,yout: 输出场的坐标
   
            % Define wavevector coordinates in frequency domain
          
            [Ny0,Nx0] = size(Uin);
            Ny=5*Ny0; Nx=5*Nx0;
            Uin_pad=zeros(Ny,Nx);
            Uin_pad(Ny/2-Ny0/2:Ny/2+Ny0/2-1,Nx/2-Nx0/2:Nx/2+Nx0/2-1)=Uin;
            dfx = 1 / (Nx * dx);
            dfy = 1 / (Ny * dx);
            fx = (-Nx/2 : Nx/2-1) * dfx;
            fy = (-Ny/2 : Ny/2-1) * dfy;
            [FX, FY] = meshgrid(fx, fy);
            k=2*pi/lambda;
            KX=2*pi*FX; KY=2*pi*FY;

            % Define transfer function
            KZ=sqrt(k^2-KX.^2 -KY.^2);
            H = exp(1j*KZ*z);
         

            % Compute Fourier transform of input field
            Uin_ft = fftshift(fft2(Uin_pad));
%             figure,imshow(abs(Uin_ft),[]);
            
            % Propagate input field in frequency domain
            Uout_ft = H .* Uin_ft;
            
            % Compute inverse Fourier transform to obtain output field
            Uout = ifft2(ifftshift(Uout_ft));
            Uout = Uout(Ny/2-Ny0/2:Ny/2+Ny0/2-1,Nx/2-Nx0/2:Nx/2+Nx0/2-1);
            % Define output coordinates
            xout = dx*(-Nx0/2:Nx0/2-1);
            yout = dx*(-Ny0/2:Ny0/2-1);
%             figure,imshow(abs(Uout),[]);

         end
    
          function [Uout,xout,yout]=prop_lens(Uin,lambda,f,dx)
            [Ny,Nx]=size(Uin);
            k=2*pi/lambda;
            x=dx*(-Nx/2:Nx/2-1);
            y=dx*(-Ny/2:Ny/2-1);
            [X,Y]=meshgrid(x,y);
            D=Nx*dx;
            t=(X.^2+Y.^2<=(D/2)^2).*exp(-1j*k/(2*f)*(X.^2+Y.^2));
            [Uout,xout,yout]=OpticUtil.prop_Angular_Spectrum(Uin.*t,lambda,f,dx);
          end

          function phase=GS(image_in,lambda,z,slm_h,slm_p,iter_num,verbose)
            arguments
                image_in
                lambda
                z
                slm_h
                slm_p
                iter_num = 100
                verbose = 0
            end
            img=im2double(image_in);
            % GS phase retrieval for lens imaging
            [Nx,Ny]=size(img);
            x=slm_p*(-Nx/2:Nx/2-1);
            y=slm_p*(-Ny/2:Ny/2-1);
            [X,Y]=meshgrid(x,y);
            A0 = (X.^2 + Y.^2 <= (slm_h/2).^2); % 光束直径限制
            a = sum(sum(A0.^2))/sum(sum(image_in.^2))/(lambda*z).^2; % 能量守恒
        
            img = image_in*sqrt(a);
            phase=rand(size(A0))*2*pi; % 初始随机相位
            for i = 1:iter_num
                Af = A0;
                f0 = Af.*exp(1i.*phase); % 振幅置1保留相位
                g0 = fftshift(fft2(f0)); 
                ang0 = angle(g0); % 取出相位与振幅
                Ampg = img;
                g1=Ampg.*exp(1i.*ang0); % 目标振幅替换，保留相位
                f1=ifft2(fftshift(g1));    
                phase=angle(f1); % Ampf = abs(f1);%取出一次迭代后相位
            end
            phase = phase + pi; % angle返回-pi~pi，转换到0~2pi
            phase = phase.*A0;
            if verbose
                figure,imshow(phase,[]);title('phase extracted');
            end
          end       


        function img_expand=expand_img(img,factor)
            factor=round(factor);
            img_expand=zeros(factor*size(img));
            for i=0:factor-1
                for j=0:factor-1
                    img_expand(1+i:factor:size(img_expand,1)+i,1+j:factor:size(img_expand,2)+j)=img;
                end
            end
        end

        function img_pad=pad_img(img,pad_sz)
            img_pad=zeros(pad_sz);
            img_pad(pad_sz(1)/2-size(img,1)/2+1:pad_sz(1)/2+size(img,1)/2,...
                    pad_sz(2)/2-size(img,2)/2+1:pad_sz(2)/2+size(img,2)/2)=img;
        end
        
       function amp_image=lee_hologram(phase_in,alpha)
            
            dim_in = size(phase_in);
            amp_image = zeros(dim_in);
            y=1:dim_in(1);
            x=1:dim_in(2);
            [X,Y]=meshgrid(x,y);
            x_y=X-Y;

%              alpha =0.85;%carrier frequency must be big enough to separate -1st order from 0th order
             amp_holo = 0.5*(1+cos(2*pi*(x_y)*alpha-phase_in)); % amplitude hologram
             amp_image(amp_holo>0.5)=1;% binary amplitude hologram
              
        end
    end
end
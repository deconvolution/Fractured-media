%% Parameters
clear all;
close all;
load('kl.mat');
load('C_b.mat');
%% Dimensions
dims.dh=6; % Spacial grid step
dims.dt=10^-2; % [s]
dims.nt=150; % Amount of time steps
dims.ns=150;
b=zeros(4,1);

b(1)=40;
b(2)=b(1)+fix((kl.i(2)-kl.i(1))/dims.dh);
b(3)=fix((kl.i(3)-kl.i(2))/dims.dh)+b(2);
b(4)=fix((kl.i(4)-kl.i(3))/dims.dh)+b(3);

dims.nz=b(4)+50;
dims.nx=2*dims.nz;
%% elastic layer
vv0=kl.vp0;
v0=kl.vp0.*sqrt(1+2*kl.del);
eta0=(kl.eps-kl.del)./(1+2*kl.del);

eta=zeros(dims.nz,dims.nx);
v=eta;
vv=eta;

% Block 1
eta(1:b(2),:)=eta0(1);
v(1:b(2),:)=v0(1);
vv(1:b(2),:)=vv0(1);
% Block 2
eta(b(2):b(3)+5,:)=eta0(2);
v(b(2):b(3)+5,:)=v0(2);
vv(b(2):b(3)+5,:)=vv0(2);
% Block 4
eta(b(4)-5:end,:)=eta0(4);
v(b(4)-5:end,:)=v0(4);
vv(b(4)-5:end,:)=vv0(4);
%% cracked layer
crack_region=[b(3),b(4)];
phi=zeros(dims.nz,dims.nx);
rho=phi;
rhof=phi;
alpf=phi;
bet=phi;

phi(:)=.0001;
rho(:)=kl.rho(3);
vv(b(3):b(4),:)=kl.vp0(3);
rhof(:)=600;
alpf(:)=1200;
bet(:)=kl.vs0(3);
%% Model dimensions
vb=ones(size(rho))*2500;
lp=40;
d0=ones(size(rho))*-3*2700*log(.00001)/2/dims.dh/40;
dims.mz=lp+1:dims.nz-lp-2;
dims.mx=lp+3:dims.nx-lp-2;
%% Source and source signals
dims.sx=fix(dims.nx)/2;
dims.sz=b(1)+1;
sn=length(dims.sx);
singles=rickerWave(10,dims);
plot(singles);
%% Receiver locations
dims.rx=min(dims.mx):max(dims.mx);
dims.rz=min(dims.mz)*ones(1,length(dims.rx));
%% source
fs=1/dims.dt;
L=dims.nt;
n=L;
f=fs*(0:(n/2))/n;
s=zeros([dims.nz*dims.nx,1]);
source_freq=fft(singles,n)/(n/2);
source_freq2=source_freq(1:n/2+1);
plot(f,abs(source_freq2));
%% check source
source_time=ifft(source_freq2,n,1)*n;
figure(3)
subplot(2,1,1)
plot(dims.dt:dims.dt:dims.dt*L,real(source_time));
subplot(2,1,2)
plot(dims.dt:dims.dt:dims.dt*L,real(singles));
%%
s_diff=diff(abs(source_freq2));
s_diff=[s_diff(1);s_diff];
s_lim=find(abs(source_freq2)<.01*max(abs(source_freq2)) & s_diff<0);

s_lim2=s_lim(1);
f_range=1:s_lim2;
f2=f(f_range);
ome=2*pi*f2;
s=zeros(dims.nx*dims.nz,length(ome));
s((dims.sx-1)*dims.nz+dims.sz,:)=source_freq2(f_range);
%%
[p,AA]=Sch_solver(vv,v,eta,phi,rho,rhof,alpf,bet,vb,d0,lp,dims,ome,s,crack_region);
%% frequency

for i=1:size(p,3)
    figure(1)
    
    subplot(2,1,1)
    imagesc(real(p(:,:,i)));
    %grid on;
    colorbar;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    title({[ '\omega=' num2str(ome(i)) 'rad/s' ],['p']});
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    axis on;
    legend([ax2,ax,ax3,ax4],'source','receiver','PML boundary','block boundary','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    
    shg;
    pause(.1);
    
end
pt=ifft(p,dims.nt,3)*dims.nt;
%% time
rec2=pt(min(dims.mz),dims.mx,:);
rec=reshape(rec2,[size(rec2,2),size(rec2,3)])';
%%

for i=1:10:size(pt,3)
    
    figure(2)
    %set(gcf,'Visible','on');
    %set(gcf,'position',[80,80,1500,600]);
    subplot(2,1,1)
    imshow(real(pt(:,:,i)),.5*[min(real(pt(:))),max(real(pt(:)))]);
    colorbar;
    axis on;
    hold on;
    ax=plot([min(dims.rx),max(dims.rx)],[dims.rz(1),dims.rz(1)],'color','red','linewidth',2);
    hold on;
    ax2=scatter(dims.sx,dims.sz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([dims.nx-lp-1,dims.nx-lp-1],[lp+1,dims.nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[lp+1,lp+1],'color','blue');
    hold on;
    ax3=plot([lp+1,dims.nx-lp-1],[dims.nz-lp-1,dims.nz-lp-1],'color','blue');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(2),b(2)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(3),b(3)],'color','yellow');
    hold on;
    ax4=plot([min(dims.mx),max(dims.mx)],[b(4),b(4)],'color','yellow');
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['t=' num2str(i*dims.dt) 's'],['p']});
    legend([ax2,ax,ax3,ax4],'source','receiver','PML boundary','block boundary','Location',[0.5,0.03,0.005,0.005],'orientation','vertical');
    shg;
    pause(.1);
    
end

%grid on;
%% Model 2 oil
phi(:)=.0001;
rhof(:)=600;
alpf(:)=1200;
[p,AA]=Sch_solver(vv,v,eta,phi,rho,rhof,alpf,bet,vb,d0,lp,dims,ome,s,crack_region);
pt=ifft(p,dims.nt,3);
rec2=pt(min(dims.mz),dims.mx,:);
rec=reshape(rec2,[size(rec2,2),size(rec2,3)])';
%save('C:\Users\Yi\OneDrive\master thesis\Schoenberg\m2_oil_pt','pt','-v7.3');
%save('C:\Users\Yi\OneDrive\master thesis\Schoenberg\m2_oil_p','p','-v7.3');
%% Model 2 water
phi(:)=.0001;
rhof(:)=1000;
alpf(:)=1500;
[p,AA]=Sch_solver(vv,v,eta,phi,rho,rhof,alpf,bet,vb,d0,lp,dims,ome,s,crack_region);
pt=ifft(p,dims.nt,3)*dims.nt;
rec2=pt(min(dims.mz),dims.mx,:);
rec=reshape(rec2,[size(rec2,2),size(rec2,3)])';
%save('C:\Users\Yi\OneDrive\master thesis\Schoenberg\m2_water_pt','pt','-v7.3');
%save('C:\Users\Yi\OneDrive\master thesis\Schoenberg\m2_water_p','p','-v7.3');
%% Model 2 gas
phi(:)=.0001;
rhof(:)=600;
alpf(:)=400;
[p,AA]=Sch_solver(vv,v,eta,phi,rho,rhof,alpf,bet,vb,d0,lp,dims,ome,s,crack_region);
pt=ifft(p,dims.nt,3);
rec2=pt(min(dims.mz),dims.mx,:);
rec=reshape(rec2,[size(rec2,2),size(rec2,3)])';
%save('C:\Users\Yi\OneDrive\master thesis\Schoenberg\m2_gas_pt','pt','-v7.3');
%save('C:\Users\Yi\OneDrive\master thesis\Schoenberg\m2_gas_pt','p','-v7.3');
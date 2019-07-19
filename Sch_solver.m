function [p,AA]=Sch_solver(vv,v,eta,phi,rho,rhof,alpf,bet,vb,d0,lp,dims,ome,s,crack_region)
%% elastic layer
AA=sparse(dims.nz*dims.nx,dims.nz*dims.nx);
p=zeros(dims.nz,dims.nx,length(ome));
D1=zeros(dims.nz,dims.nx);
D2=D1;
D3=D1;
%% porous layer
B1=zeros(dims.nz,dims.nx);
B2=B1;
B3=B2;
B4=B3;
B5=B4;
alppl=2*(1-bet.^2./vv.^2).^.5.*bet;
hs=1-phi;
%%
tic;
for l=2:length(ome)
    %% assign d
    dx=zeros(dims.nz,dims.nx);
    dz=zeros(dims.nz,dims.nx);
    % top
    for i=3:lp+2
        for j=lp+3:dims.nx-lp-2
            dx(i,j)=0;
            dz(i,j)=d0(i,j)*((lp+3-i)/lp)^2;
        end
    end
    % bottom
    for i=dims.nz-lp-1:dims.nz-2
        for j=lp+3:dims.nx-lp-2
            dx(i,j)=0;
            dz(i,j)=d0(i,j)*((i+2-(dims.nz-lp))/lp)^2;
        end
    end
    % left
    for i=lp+3:dims.nz-lp-2
        for j=3:lp+2
            dz(i,j)=0;
            dx(i,j)=d0(i,j)*((lp+3-j)/lp)^2;
        end
    end
    % right
    for i=lp+2:dims.nz-lp-2
        for j=dims.nx-lp-1:dims.nx-2
            dz(i,j)=0;
            dx(i,j)=d0(i,j)*((j+2-(dims.nx-lp))/lp)^2;
        end
    end
    % upper left
    for i=3:lp+2
        for j=3:lp+2
            dz(i,j)=d0(i,j)*((lp+3-i)/lp)^2;
            dx(i,j)=d0(i,j)*((lp+3-j)/lp)^2;
        end
    end
    % upper right
    for i=3:lp+2
        for j=dims.nx-lp-1:dims.nx-2
            dz(i,j)=d0(i,j)*((lp+3-i)/lp)^2;
            dx(i,j)=d0(i,j)*((j+2-(dims.nx-lp))/lp)^2;
        end
    end
    % lower left
    for i=dims.nz-lp-1:dims.nz-2
        for j=3:lp+2
            dz(i,j)=d0(i,j)*((i+2-(dims.nz-lp))/lp)^2;
            dx(i,j)=d0(i,j)*((lp+3-j)/lp)^2;
        end
    end
    % lower right
    for i=dims.nz-lp-1:dims.nz-2
        for j=dims.nx-lp-1:dims.nx-2
            dz(i,j)=d0(i,j)*((i+2-(dims.nz-lp))/lp)^2;
            dx(i,j)=d0(i,j)*((j+2-(dims.nx-lp))/lp)^2;
        end
    end
    gx=1+dx/1i./ome(l);
    gz=1+dz/1i./ome(l);
    gx_x=diff(gx,1,2)/dims.dh;
    gz_z=diff(gz,1,1)/dims.dh;
    gx_x=[-gx_x(:,1:lp+2),zeros([size(gx_x,1),1]),gx_x(:,lp+3:end)];
    gx_x(:,[1,2])=0;
    gx_x(:,[end,end-1])=0;
    gz_z=[-gz_z(1:lp+2,:);zeros([1,size(gz_z,2)]);gz_z(lp+3:end,:)];
    gz_z([1,2],:)=0;
    gz_z([end,end-1],:)=0;
    %% B
    B1=ome(l)^4./alppl.^2.*(phi./rhof./alpf.^2+(1-phi)./rho./vv.^2);
    B2=(ome(l)^2./alppl.^2.*(phi./rhof+(1-phi)./rho)+ome(l)^2*(phi./rhof./alpf.^2+(1-phi)./rho./vv.^2)-ome(l)^2.*hs./vv.^2./rho)./gx.^2;
    B3=ome(l)^2./alppl.^2./(phi.*rhof+(1-phi).*rho)./gz.^2;
    B4=((phi./rhof+(1-phi)./rho)-hs./rho)./gx.^4;
    B5=1./(phi.*rhof+(1-phi).*rho)./gx.^2./gz.^2;
    %% D
    D1=(1+2*eta).*v.^2./gx.^2;
    D2=vv.^2./gz.^2;
    D3=2*eta.*v.^2.*vv.^2/ome(l)^2.*gx.^2.*gz.^2;
    %% domain interior
    for i=3:dims.nz-2
        for j=3:dims.nx-2
            if i<crack_region(1) || i>crack_region(2)
                % i,j
                AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2*D1(i,j)/dims.dh^2-2*D2(i,j)/dims.dh^2+4*D3(i,j)/dims.dh^4+ome(l)^2;
                % i,j+1
                AA((j-1)*dims.nz+i,(j)*dims.nz+i)=D1(i,j)/dims.dh^2-2*D3(i,j)/dims.dh^4;
                % i,j-1
                AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=D1(i,j)/dims.dh^2-2*D3(i,j)/dims.dh^4;
                % i+1,j
                AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=D2(i,j)/dims.dh^2-2*D3(i,j)/dims.dh^4;
                % i-1,j
                AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=D2(i,j)/dims.dh^2-2*D3(i,j)/dims.dh^4;
                % i+1,j+1
                AA((j-1)*dims.nz+i,(j)*dims.nz+i+1)=D3(i,j)/dims.dh^4;
                % i-1,j+1
                AA((j-1)*dims.nz+i,(j)*dims.nz+i-1)=D3(i,j)/dims.dh^4;
                % i+1,j-1
                AA((j-1)*dims.nz+i,(j-2)*dims.nz+i+1)=D3(i,j)/dims.dh^4;
                % i-1,j-1
                AA((j-1)*dims.nz+i,(j-2)*dims.nz+i-1)=D3(i,j)/dims.dh^4;
            else
                % i,j
                AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=B1(i,j)-2*B2(i,j)/dims.dh^2-2*B3(i,j)/dims.dh^2+6*B4(i,j)/dims.dh^4+4*B5(i,j)/dims.dh^4;
                
                % i,j+1
                AA((j-1)*dims.nz+i,(j)*dims.nz+i)=B2(i,j)/dims.dh^2-4*B4(i,j)/dims.dh^4-2*B5(i,j)/dims.dh^4;
                % i,j-1
                AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=B2(i,j)/dims.dh^2-4*B4(i,j)/dims.dh^4-2*B5(i,j)/dims.dh^4;
                
                % i,j+2
                AA((j-1)*dims.nz+i,(j+1)*dims.nz+i)=B4(i,j)/dims.dh^4;
                % i,j-2
                AA((j-1)*dims.nz+i,(j-3)*dims.nz+i)=B4(i,j)/dims.dh^4;
                
                % i+1,j+1
                AA((j-1)*dims.nz+i,(j)*dims.nz+i+1)=B5(i,j)/dims.dh^4;
                % i-1,j-1
                AA((j-1)*dims.nz+i,(j-2)*dims.nz+i-1)=B5(i,j)/dims.dh^4;
                
                % i+1,j
                AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=B3(i,j)/dims.dh^2-2*B5(i,j)/dims.dh^4;
                % i-1,j
                AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=B3(i,j)/dims.dh^2-2*B5(i,j)/dims.dh^4;
                
                % i-1,j+1
                AA((j-1)*dims.nz+i,(j)*dims.nz+i-1)=B5(i,j)/dims.dh^4;
                % i+1,j-1
                AA((j-1)*dims.nz+i,(j-2)*dims.nz+i+1)=B5(i,j)/dims.dh^4;
            end
        end
    end
    %% first boundary line
    %% B2
    i=1;
    for j=2:dims.nx-1
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh^2;
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B3
    i=dims.nz;
    for j=2:dims.nx-1
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B4
    j=dims.nx;
    for i=2:dims.nz-1
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B5
    j=1;
    for i=2:dims.nz-1
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B6
    j=dims.nx;
    i=1;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% B7
    j=1;
    i=1;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% B8
    j=dims.nx;
    i=dims.nz;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% B9
    j=1;
    i=dims.nz;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% second boundary line
    %% B2
    i=2;
    for j=3:dims.nx-2
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh^2;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh^2;
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B3
    i=dims.nz-1;
    for j=3:dims.nx-2
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B4
    j=dims.nx-1;
    for i=3:dims.nz-2
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B5
    j=2;
    for i=3:dims.nz-2
        AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-1/dims.dh-1i*ome(l)/(vb(i,j));
        AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
        %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    end
    %% B6
    j=dims.nx-1;
    i=2;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% B7
    j=2;
    i=2;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% B8
    j=dims.nx-1;
    i=dims.nz-1;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    %AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %% B9
    j=2;
    i=dims.nz-1;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i)=-2/dims.dh-2i*ome(l)/(vb(i,j));
    AA((j-1)*dims.nz+i,(j)*dims.nz+i)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-2)*dims.nz+i)=1/dims.dh;
    AA((j-1)*dims.nz+i,(j-1)*dims.nz+i-1)=1/dims.dh;
    %AA((j-1)*dims.nz+i,(j-1)*dims.nz+i+1)=1/dims.dh;
    %%
    u0=AA\-s(:,l);
    for i=1:length(u0)
        p(i+dims.nz*dims.nx*(l-1))=u0(i);
    end
    fprintf('l=%f\n',l);
    toc;
end
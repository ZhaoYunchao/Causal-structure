close all; clear all;

%% Parameters
ut=0.0637541;% based on ub=1, do not modify
ub=0.8982822;% the current data is normalized by uc
Reb=3200*ub;
Ret=Reb*ut;
nu=1/Reb;

nx=256;
nz=257;
lx=4*pi;
lz=2*pi;
dx=lx/nx;
dz=lz/nz;
nt=8000;
x=(1:nx)*dx;
z=(1:nz)*dz;
dt=0.085741;% based on ub=1
dtplus=dt*ut*Ret;
timedelay_length=101;

%% Reading \tau_w
tauw=zeros(nx,nz,nt);
for l=1:nt
    filename = ['/home/zhaoyunchao/OpenChannel/twolayers_u/u_yplus20_',num2str(l,'%05d'),'.u'];
    fileID = fopen(filename);
    tauw(:,:,l)=fread(fileID,[nz,nx],'single')'/ub;
    fclose(fileID);
end
tauw=tauw/ut^2;

tauw_mean=mean(tauw,'all');
tauw=tauw-tauw_mean;

dtauwdt=(tauw(:,:,2:end)-tauw(:,:,1:end-1))/dt;
dtauwdt_mean=mean(dtauwdt,'all');
dtauwdt=dtauwdt-dtauwdt_mean;

%% Reading u and computing LKIF
u_3d_T21=zeros(63,33,timedelay_length,63);
u_3d_tau21=zeros(63,33,timedelay_length,63);
u_3d_C12=zeros(63,33,timedelay_length,63);

for j=3:65 % the 35th layer ==>  y^+=60
    u=zeros(nx,nz,nt);
    disp(['Reading j=',num2str(j)]);
    for l=1:nt
        file_u = ['E:\Jimenez new channel\Openchannel\U\re180caus.',num2str(l,'%05d'),'.uvw.u'];
        [para,yvec,field_u] = read_flow(file_u);
        u(:,:,l)=(squeeze(field_u(:,:,j)))'/ub;
    end
    disp(['Reading j=',num2str(j),' done']);
    u=u/ut;

    u2_mean=mean(u,'all');
    u=u-u2_mean;

    T21=zeros(63,33,timedelay_length);
    tau21=zeros(63,33,timedelay_length);
    C12=zeros(63,33,timedelay_length);

    for ii=1:7
        for kk=1:15
            u1=squeeze(tauw(32*ii,1+16*kk,1:end-1));
            u2=u((-31:31)+32*ii,(-15:17)+16*kk,1:end-1);
            du1dt=squeeze(dtauwdt(32*ii,1+16*kk,:));
            du2dt=dudt((-31:31)+32*ii,(-15:17)+16*kk,1:end-1);

            [T21_temp,tau21_temp,C12_temp]=Liang(u1,u2,du1dt);
            T21=T21+T21_temp;
            tau21=tau21+tau21_temp;
            C12=C12+C12_temp;
        end
        disp(['y No.',num2str(j),', x No.',num2str(ii)]);
    end
    T21=T21/105;
    tau21=tau21/105;
    C12=C12/105;

    u_3d_T21(:,:,:,j)=T21;
    u_3d_tau21(:,:,:,j)=tau21;
    u_3d_C12(:,:,:,j)=C12;
end

delta_x=(-31:31)*dx;
delta_z=(-16:16)*dz;
delta_t=(-60:40)*dt;
delta_tplus=(-60:40)*dtplus;

%% Subfunction
function [T21,tau21,C12]=Liang(u1_input,u2_input,du1dt_input)
    timedelay_length=101;
    delta_t=-60:40;
    dt=0.085741;
    
    T21=zeros(63,33,timedelay_length);
    tau21=zeros(63,33,timedelay_length);
    C12=zeros(63,33,timedelay_length);
    for m=1:timedelay_length
        if delta_t(m)<=0
            u1=repmat(shiftdim(u1_input(1-delta_t(m):end)',-1),63,33);
            u2=u2_input(:,:,1:end+delta_t(m));
            du1dt=repmat(shiftdim(du1dt_input(1-delta_t(m):end)',-1),63,33);
        else
            u1=repmat(shiftdim(u1_input(1:end-delta_t(m))',-1),63,33);
            u2=u2_input(:,:,delta_t(m)+1:end);
            du1dt=repmat(shiftdim(du1dt_input(1:end-delta_t(m))',-1),63,33);
        end
        
        C11=mean(u1.^2,3);
        C12_temp=mean(u1.*u2,3);
        C22=mean(u2.^2,3);
        C1d1=mean(u1.*du1dt,3);
        C2d1=mean(u2.*du1dt,3);
        Cd1d1=mean(du1dt.^2,3);
        
        p1=(C22.*C1d1-C12_temp.*C2d1)./(C11.*C22-C12_temp.^2);
        q1=(C11.*C2d1-C12_temp.*C1d1)./(C11.*C22-C12_temp.^2);
        noise1=dt/2*(Cd1d1+p1.^2.*C11+q1.^2.*C22-2*p1.*C1d1-2*q1.*C2d1+2*p1.*q1.*C12_temp)./C11;
        T21(:,:,m)=T21(:,:,m)+C12_temp.*q1./C11;
        tau21(:,:,m)=tau21(:,:,m)+C12_temp.*q1./C11./(abs(C12_temp.*q1./C11)+abs(p1)+abs(noise1));
        C12(:,:,m)=C12_temp./sqrt(C11.*C22);
    end
end
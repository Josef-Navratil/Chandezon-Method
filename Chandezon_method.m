clear all; clc;
%% Prepare constants
tic                     %start counting computing time
[Parameters]=setParameters();
VecWaveLength=linspace(Parameters(1).Param(1),Parameters(1).Param(2),Parameters(1).Param(3));
RefVec=zeros(Parameters(1).Param(3),1);
TranVec=zeros(Parameters(1).Param(3),1);
WaveLengthToPlot=Parameters(1).Param(11);
%% Set constants
d=Parameters(1).Param(4);                %depth of the grating            
thI=Parameters(1).Param(5)*pi/180;          %incident angle
n1=Parameters(1).Param(6);                   %refraction index of incident medium
n2=Parameters(1).Param(7);              %refraction index of transmission medium 
K=2*pi/d;                  
mu0=4*pi*1e-7;          %vaccum permeability
eps0=8.85*1e-12;        %vacuum permittivity
nTr=Parameters(1).Param(8);   %Truncation number
nDim=2*nTr+1;              %Total number of modes
tol=Parameters(1).Param(9);             %error tolerance
h=Parameters(1).Param(10);

%Z0=sqrt(mu0/eps0);     %Z0
%% Set polarization
if strcmp(Parameters(1).Polar,'TM') %TM polarization
mu=1;                   %relative permeability
eps1=n1*n1/(mu^2);      %relative permittivity of incident medium
eps2=n2*n2/(mu^2);      %relative permittivity of transmission medium
end

if strcmp(Parameters(1).Polar,'TE') %TE polarization
eps1=1;                %relative permittivity of incident medium
eps2=1;                %relative permittivity of transmitted medium
mu1=n1;                 %relative permeability of incident medium
mu2=n2;              %relative permeability of transmitted medium
end

for i=1:length(VecWaveLength); %start computing diffraction efficiencies for given wavelengths
Lam=VecWaveLength(i);

k0=2*pi/Lam;            %wavenumber

%n2=1.5;
%% Grating profile and its derivative

a_fun=@(x) k0.*(Parameters(1).Prof(x./k0));
a_diff_fun=@(x) Parameters(2).Prof(x./k0);

%a_fun=@(x) (h/2*k0)*cos((K/k0).*x);
%a_diff_fun=@(x) -(h/2*K).*sin((K/k0).*x);

cut_small=1;            %cut small array elements

%% Auxiliary constants
%m1=-28;                  %lower truncation order
%m2=26;                   %upper truncation order
alfa0=n1*k0*sin(thI);
m1=-floor(alfa0/K)-(nDim-1)/2;     %lower order for adaptive truncation
m2=-floor(alfa0/K)+(nDim-1)/2;     %upper order for adaptive truncation
%% Prepare fields
nDim=m2-m1+1;           %number of modes
alfa0=n1*k0*sin(thI);   %alpha0 from (3)
A=(alfa0*ones(1,m2-m1+1)+K.*linspace(m1,m2,nDim))./k0; %alpha field from (2)

B1=n1^2-A.^2;                       %beta_1^2 field from (4)
B2=n2^2-A.^2;                       %beta_2^2 field from (4)

if (min(abs(B1)<tol))               %When some element of B1 is equal to zero, break
    disp('System resonance, change the incident angle slightly');
    return;
end
if (min(abs(B2)<tol))               %When some element of B2 is equal to zero, break
    disp('System resonance, change the incident angle slightly')
    return;
end
SB1=sqrt(B1);                       %Eigenvalues of the original problem, incident medium
SB1_idx=(abs(imag(SB1))==0)&(real(SB1)>0); %indices of positive real propagation orders, incident medium
SB1_ind=(m1:m2);                    %indices of all modes
real_Ray1_idx=SB1_ind(SB1_idx);     %indices of real propagation numbers, positive due to radiation conditions

SB2=-sqrt(B2);                      %Eigenvalues of the original problem, incident medium
SB2_idx2=(abs(imag(SB2))==0)&(real(SB2)<0); %indices of negative real propagation orders, transmission medium
SB2_ind2=(m1:m2);                   %indices of all modes
real_Ray2_idx=SB2_ind2(SB2_idx2);   %indices of all modes
%% Evaluate the FFT of function a_diff
a_diff_vec=F_series_gen(a_diff_fun,12,k0*d,nDim);
a_mat=toeplitz(a_diff_vec);


%% 
[V1,rho1,V2,rho2]=EigChand(A,B1,B2,a_mat,nDim);
%clear ChandM1 ChandM2 BB1 BB2

%% Sort eigenvalues
[real_eig1p,real_eig2n,imag_eig1p,imag_eig2n,imag_Vec1p,imag_Vec2n]=SortEigenvaluesChand(rho1,V1,rho2,V2,tol,nDim);
%clear rho1 rho2 V1 V2      %If the program
%has problems with memory, delete superfluous variables

%% Assemble F-matrices
b0=sqrt(B1(1-m1));

[F_in,FRN,FRP]=GenerateFFieldsChand(a_fun,b0,nDim,k0,d,m1,m2,real_Ray2_idx,real_Ray1_idx,SB1,SB2);

%% Assemble G matrices
[G_RP,G_RN,G_in,G_P,G_N]=GenerateGFieldsChand(b0,a_mat,real_eig1p,real_eig2n,SB1,SB2,real_Ray1_idx,real_Ray2_idx,m1,m2,nDim,imag_eig1p,imag_eig2n,FRP,FRN,F_in,imag_Vec1p,imag_Vec2n,A,eps1,eps2);

%% Assemble matrix for matching boundary conditions and solve the linear system
MatBC=[FRP imag_Vec1p -FRN -imag_Vec2n;G_RP G_P -G_RN -G_N];
VecBC=[-F_in; -G_in];
RVec=MatBC\VecBC;


%% Print Efficiencies
etaR=zeros(1,length(real_eig1p));
etaT=zeros(1,length(real_eig2n));
for k=min(real_Ray1_idx):max(real_Ray1_idx)
    etaR(k-min(real_Ray1_idx)+1)=(sqrt(B1(k-m1+1))/b0)*(abs(RVec(k-min(real_Ray1_idx)+1)))^2;
end
for k=min(real_Ray2_idx):max(real_Ray2_idx)
    etaT(k-min(real_Ray2_idx)+1)=(eps1/eps2)*(sqrt(B2(k-m1+1))/b0)*(abs(RVec(k-min(real_Ray2_idx)+1+nDim)))^2;
end
if((VecWaveLength(i)<=WaveLengthToPlot)&&(WaveLengthToPlot<=VecWaveLength(i+1)))
    disp('Plotting intensity distribution');
    %% Plot fields
    nPoints=400;
    xx=linspace(0,d*k0,nPoints);
    yy=linspace(-0.7*h*k0,h*k0,nPoints);
    [X,Y]=meshgrid(xx,yy);
    clear yy

    FIn=exp(1i*alfa0/k0*X-1i*b0*Y).*(a_fun(X)<=Y);
    FRPlus=zeros(nPoints);          %preallocate fields for positive propagation orders in superstrate medium
    FRNeg=zeros(nPoints);           %preallocate fields for negative propagation orders in substrate medium
    FRPlusIm=zeros(nPoints);        %preallocate fields for positive evanescent orders in superstrate medium
    FRNegIm=zeros(nPoints);         %preallocate fields for negative evanescent orders in substrate medium

    for mm=1:length(real_eig1p)
        FRPlus=(FRPlus + RVec(mm).*exp(1i*A(mm+min(real_Ray1_idx)-m1).*X + 1i*SB1(mm+min(real_Ray1_idx)-m1).*Y)).*(a_fun(X)<=Y);
    end
 
    for mm=m1:m2
        for kk=1:length(imag_eig1p)
            FRPlusIm=FRPlusIm + (exp(1i*A(mm-m1+1).*X).*((RVec(kk+length(real_eig1p))*imag_Vec1p(mm-m1+1,kk)).*exp(1i*imag_eig1p(kk).*(Y-a_fun(X))))).*(a_fun(X)<=Y);
        end
    end
    disp('Keep waiting, we are working on it')
    for mm=1:length(real_eig2n)
       FRNeg=FRNeg + (RVec(mm+length(real_eig1p)+length(imag_eig1p)).*exp(1i*A(mm+min(real_Ray2_idx)-m1).*X + 1i*SB2(mm+min(real_Ray2_idx)-m1).*Y)).*(a_fun(X)>Y);
    end

    for mm=m1:m2
        for kk=1:length(imag_eig2n)
           FRNegIm=FRNegIm + (exp(1i*A(mm-m1+1).*X).*((RVec(kk+length(real_eig1p)+length(imag_eig1p)+length(real_eig2n))*imag_Vec2n(mm-m1+1,kk)).*exp(1i*imag_eig2n(kk).*(Y-a_fun(X))))).*(a_fun(X)>Y);
        end
    end
    disp('It is almost done')
    Z=double(abs(FIn+FRPlus+FRPlusIm+FRNeg+FRNegIm));
    %Z=double(real(FIn+FRPlus+FRPlusIm+FRNeg+FRNegIm));
    pcolor(X,Y,Z);
    shading flat;
    colorbar;
    colormap jet;
    hold on
    fy=a_fun(xx);
    plot(xx,fy,'w-','LineWidth',1.5);
    disp('Computation paused, press any key to continue')
    pause

end

clearvars -except etaR etaT Lam WaveLengthToPlot VecWaveLength Parameters d h thI n1 n2 K mu0 eps0 nTr tol alfa0 mu eps1 eps2 mu1 mu2 nDim real_Ray1_idx i RefVec
RAmp=(real_Ray1_idx==0);
RefVec(i)=etaR(RAmp);
clear real_Ray1_idx
toc
end
plot(VecWaveLength,RefVec);
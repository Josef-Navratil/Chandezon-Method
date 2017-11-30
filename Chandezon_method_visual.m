%% Implementation of Chandezon method for relief gratings with visualisation of the fields
clear all; %clc;
%% Set constants
tic                     %start counting computing time
Lam=400e-9;             %wavelength of incident light
h=600e-9;                 %depth of the grating
d=300e-9;               %period of the grating
K=2*pi/d;               
mu0=4*pi*1e-7;          %vaccum permeability
eps0=8.85*1e-12;        %vacuum permittivity
thI=15*pi/180;          %incident angle
n1=1;                   %refraction index of incident medium
n2=1+5*1i;              %refraction index of transmission medium 
k0=2*pi/Lam;            %wavenumber
%n2=1.5;
%% Grating profile and its derivative

%Asymetrical profile
%a_fun=@(x) (0.1*Lam*k0)*cos((K/k0).*x) + (0.02*Lam*k0)*cos((2*K/k0).*x - 5*pi/9);
%a_diff_fun=@(x) -(0.1*Lam*K).*sin((K/k0).*x) - (0.02*Lam*2*K).*(sin((2*K/k0).*x - 5*pi/9));

%Symmetrical profile
a_fun=@(x) (h/2*k0)*cos((K/k0).*x);
a_diff_fun=@(x) -(h/2*K).*sin((K/k0).*x);

%Trivial profile
%a_fun=@(x) 0.*x;
%a_diff_fun=@(x) 0.*x;
cut_small=1;            %cut small array elements
%% TM Polarization
mu=1;                   %relative permeability
eps1=n1*n1/(mu^2);      %relative permittivity of incident medium
eps2=n2*n2/(mu^2);      %relative permittivity of transmission medium
%% TE Polarization
%eps1=1;                %relative permittivity of incident medium
%eps2=1;                %relative permittivity of transmitted medium
%mu1=1;                 %relative permeability of incident medium
%mu2=1.5;              %relative permeability of transmitted medium
%% Auxiliary constants
m1=-25;                  %lower truncation order
m2=25;                   %upper truncation order
tol=10e-10;             %error tolerance
%Z0=sqrt(mu0/eps0);     %Z0
%m1=-floor(alfa0/K)-(nDim-1)/2;     %lower order for adaptive truncation
%m2=-floor(alfa0/K)+(nDim-1)/2;     %upper order for adaptive truncation
%Here ends the input parameters - don't touch anything below :-)
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


%% Assemble matrix and compute eigenvalues
IB1=diag(B1)\eye(nDim);
IB2=diag(B2)\eye(nDim);
AUX_mat=(eye(nDim) + a_mat^2);

ChandM1=[-IB1*(diag(A)*a_mat + a_mat*diag(A))  (IB1*AUX_mat); eye(nDim),zeros(nDim)];     %eigenvalue matrix from (12), incident medium
ChandM2=[-IB2*(diag(A)*a_mat + a_mat*diag(A))  (IB2*AUX_mat); eye(nDim),zeros(nDim)];     %eigenvalue matrix from (12), transmission medium

[V1,rho1]=eig(ChandM1);             %Eigenvalues and eigenvectors of ChandM1
[V2,rho2]=eig(ChandM2);             %Eigenvalues and eigenvectors of ChandM2

rho1=1./diag(rho1);                 %Eigenvalues of ChandM1
rho2=1./diag(rho2);                 %Eigenvalues of ChandM2
%clear ChandM1 ChandM2 BB1 BB2


%% Sort eigenvalues
real_eig1p_ind=(abs(imag(rho1))<tol)&(real(rho1)>tol);    %indices of eigenvalues with positive real part, incident medium
real_eig2n_ind=(abs(imag(rho2))<tol)&(real(rho2)<-tol);    %indices of eigenvalues with negative real part, transmission medium
real_eig1p=rho1(real_eig1p_ind);                        %positive real eigenvalues, incident medium
real_eig2n=rho2(real_eig2n_ind);                        %negative real eigenvalues, transmission medium
[sort_real1,idx1]=sort(real_eig1p,'descend');           %sort eigenvalues descendend, incident medium
[sort_real2,idx2]=sort(real_eig2n,'descend');           %sort eigenvalues descendend, transmission medium
s_real_Vec1p=V1(1:nDim,real_eig1p_ind);                 %eigenvectors of positive real eigenvalues, incident medium
s_real_Vec2n=V2(1:nDim,real_eig2n_ind);                 %eigenvectors of negative real eigenvalues, transmission medium
real_eig1p=real(real_eig1p(idx1));                      %sorted positive real eigenvalues, incident medium
real_eig2n=real(real_eig2n(idx2));                      %sorted real negative eigenvalues, transmission medium
real_Vec1p=s_real_Vec1p(:,idx1');                       %sorted real eigenvectors, incident medium
real_Vec2n=s_real_Vec2n(:,idx2');                       %sorted real eigenvectors, transmission medium
imag_eig1p_ind=(imag(rho1)>tol);                        %indices of eigenvalues with positive imaginary part, incident medium
imag_eig2n_ind=(imag(rho2)<-tol);                       %infices of eigenvalues with negative imaginary pert, transmission medium
imag_eig1p=rho1(imag_eig1p_ind);                        %eigenvalues with positive imaginary part, incident medium
imag_eig2n=rho2(imag_eig2n_ind);                        %eigenvalues with negative imaginary pert, transmission medium
[sort_imag1,idx3]=sort(abs(imag(imag_eig1p)),'ascend'); %sort eigenvalues with with positive imaginary part, incident medium
[sort_imag2,idx4]=sort(abs(imag(imag_eig2n)),'ascend'); %sort eigenvalues with negative imaginary pert, transmission medium
imag_eig1p=imag_eig1p(idx3);
imag_eig2n=imag_eig2n(idx4);
s_imag_Vec1p=V1(1:nDim,imag_eig1p_ind);
s_imag_Vec2n=V2(1:nDim,imag_eig2n_ind);
imag_Vec1p=s_imag_Vec1p(:,idx3);
imag_Vec2n=s_imag_Vec2n(:,idx4);
FPVec1=[real_Vec1p s_imag_Vec1p];
FNVec2=[real_Vec2n s_imag_Vec2n];
rho_1p=[real_eig1p; imag_eig1p];
rho_2n=[real_eig2n; imag_eig2n];
%clear rho1 rho2 rhopind rhonind rhop2ind rhon2ind V1 V2      %If the program
%has problems with memory, delete superfluous variables


%% Assemble matrix L_m
b0=sqrt(B1(1-m1));
%Generate field F_in
LP_fun=@(x) exp(-(1i*b0).*a_fun(x));
LN_fun=@(x) exp((1i*b0).*a_fun(x));
F0=F_series_gen(LP_fun,12,k0*d,nDim);
F1=F_series_gen(LN_fun,12,k0*d,nDim);
F_in=[flip(F1(2:-m1+1));conj(F0(1:m2+1))];
FRN=zeros(nDim,length(real_Ray2_idx));
FRP=zeros(nDim,length(real_Ray1_idx));
          

for M=min(real_Ray1_idx):max(real_Ray1_idx)    
    LPn_fun=@(x) exp((1i*SB1(M-m1+1)).*a_fun(x));
    fm=12;%exponent to number of fourier modes
    FRP0=F_series_gen(LPn_fun,fm,k0*d,2^fm);
    FRP(:,M-min(real_Ray1_idx)+1)=[conj(FRP0((2^fm+m1+1-M):(2^fm)));conj(FRP0(1:(m2+1-M)))];
end


for M=min(real_Ray2_idx):max(real_Ray2_idx)    
    LNn_fun=@(x) exp((1i*SB2(M-m1+1)).*a_fun(x));
    fm=12;%exponent to number of fourier modes
    FRN0=F_series_gen(LNn_fun,fm,k0*d,2^fm);
    FRN(:,M-min(real_Ray2_idx)+1)=[conj(FRN0((2^fm+m1+1-M):(2^fm)));conj(FRN0(1:(m2+1-M)))];
end


%% Assemble G matrices
G_RP=zeros(nDim,length(real_eig1p));
G_RN=zeros(nDim,length(real_eig2n));
G_in=zeros(nDim,1);
for M=1:nDim
    for N=1:length(real_eig1p)
        for NN=1:nDim
            G_RP(M,N)=G_RP(M,N)+ (a_mat(M,NN)*A(NN) - AUX_mat(M,NN)*SB1(N+min(real_Ray1_idx)-m1))*FRP(NN,N);
        end
    end
end

for M=1:nDim
    for N=1:length(real_eig2n)
        for NN=1:nDim
           G_RN(M,N)=G_RN(M,N)+ (a_mat(M,NN)*A(NN) - AUX_mat(M,NN)*SB2(N+min(real_Ray2_idx)-m1))*FRN(NN,N);
        end
    end
end

for M=1:nDim
    for NN=1:nDim
            G_in(M)=G_in(M)+ (a_mat(M,NN)*A(NN) + AUX_mat(M,NN)*b0)*F_in(NN);
    end
end

G_P=zeros(nDim,length(imag_eig1p));
G_N=zeros(nDim,length(imag_eig2n));
for M=1:nDim
    for N=1:length(imag_eig1p)
        for NN=1:nDim
            G_P(M,N)=G_P(M,N)+ (a_mat(M,NN)*A(NN) - AUX_mat(M,NN)*(imag_eig1p(N)))*imag_Vec1p(NN,N);
        end
    end
end

for M=1:nDim
    for N=1:length(imag_eig2n)
        for NN=1:nDim
            G_N(M,N)=G_N(M,N)+ (a_mat(M,NN)*A(NN) - AUX_mat(M,NN)*(imag_eig2n(N)))*imag_Vec2n(NN,N);
        end
    end
end

G_RP=(1/eps1).*G_RP;
G_RN=(1/eps2).*G_RN;
G_in=(1/eps1).*G_in;
G_P=(1/eps1).*G_P;
G_N=(1/eps2).*G_N;


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
toc
%% Plot fields
nPoints=200;
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
toc

for mm=m1:m2
    for kk=1:length(imag_eig1p)
        FRPlusIm=FRPlusIm + (exp(1i*A(mm-m1+1).*X).*((RVec(kk+length(real_eig1p))*imag_Vec1p(mm-m1+1,kk)).*exp(1i*imag_eig1p(kk).*(Y-a_fun(X))))).*(a_fun(X)<=Y);
    end
end
toc
for mm=1:length(real_eig2n)
    FRNeg=FRNeg + (RVec(mm+length(real_eig1p)+length(imag_eig1p)).*exp(1i*A(mm+min(real_Ray2_idx)-m1).*X + 1i*SB2(mm+min(real_Ray2_idx)-m1).*Y)).*(a_fun(X)>Y);
end

for mm=m1:m2
    for kk=1:length(imag_eig2n)
        FRNegIm=FRNegIm + (exp(1i*A(mm-m1+1).*X).*((RVec(kk+length(real_eig1p)+length(imag_eig1p)+length(real_eig2n))*imag_Vec2n(mm-m1+1,kk)).*exp(1i*imag_eig2n(kk).*(Y-a_fun(X))))).*(a_fun(X)>Y);
    end
end

Z=double(abs(FIn+FRPlus+FRPlusIm+FRNeg+FRNegIm));

pcolor(X,Y,Z);
shading flat;
colorbar;
colormap jet;
hold on
fy=a_fun(xx);
plot(xx,fy,'w-','LineWidth',1.5);
%{
FIn=@(x,y) exp(1i*alfa0/k0.*x-1i*b0.*y);
%FIn=FIn(X,Y);
FRPlus=@(x,y) 0;
FRPlusIm=@(x,y) 0;
FRNeg=@(x,y) 0;
FRNegIm=@(x,y) 0;

for mm=1:length(real_eig1p)
    FRPlus=@(x,y) FRPlus(x,y) + RVec(mm).*exp(1i*A(mm+min(real_Ray1_idx)-m1).*x + 1i*SB1(mm+min(real_Ray1_idx)-m1).*y);
end
%FRPlus=FRPlus(X,Y);
for mm=m1:m2
    for kk=1:length(imag_eig1p)
        FRPlusIm=@(x,y) FRPlusIm(x,y) + exp(1i*A(mm-m1+1).*x).*((RVec(kk+length(real_eig1p))*imag_Vec1p(mm-m1+1,kk)).*exp(1i*imag_eig1p(kk).*(y-a_fun(x))));
    end
end
%FRPlusIm=FRPlusIm(X,Y);
for mm=1:length(real_eig2n)
    FRNeg=@(x,y) FRNeg(x,y) + RVec(mm+length(real_eig1p)+length(imag_eig1p)).*exp(1i*A(mm+min(real_Ray2_idx)-m1).*x + 1i*SB2(mm+min(real_Ray2_idx)-m1).*y);
end
%FRNeg=FRNeg(X,Y);
for mm=m1:m2
    for kk=1:length(imag_eig2n)
        FRNegIm=@(x,y) FRNegIm(x,y) + exp(1i*A(mm-m1+1).*x).*((RVec(kk+length(real_eig1p)+length(imag_eig1p)+length(real_eig2n))*imag_Vec2n(mm-m1+1,kk)).*exp(1i*imag_eig2n(kk).*(y-a_fun(x))));
    end
end
%FRNegIm=FRNegIm(X,Y);
toc

FX=@(x,y) (FIn(x,y) + FRPlus(x,y)+FRPlusIm(x,y)).*(y>a_fun(x))+ (FRNeg(x,y)+FRNegIm(x,y)).*(y<=a_fun(x));

%FX=FRPlus+FRPlusIm+FRNeg+FRNegIm;

FFX=FX(X,Y);
Z=double(abs(FFX));
pcolor(X,Y,Z);
shading flat;
colorbar;
colormap jet;
hold on
fy=a_fun(xx);
plot(xx,fy,'w-','LineWidth',0.3);
%}
%clearvars -except etaR etaT Lambda
toc
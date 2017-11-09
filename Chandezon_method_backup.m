clear all; clc;
%% Prepare constants
tic                     %start counting computing time
[Parameters]=setParameters();
VecWaveLength=linspace(Parameters(1).Param(1),Parameters(1).Param(2),Parameters(1).Param(3));
RefVec=zeros(Parameters(1).Param(3),1);
TranVec=zeros(Parameters(1).Param(3),1);
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




for i=1:length(VecWaveLength);
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
AUX_mat=(eye(nDim) + a_mat^2);
for M=1:nDim
    for N=1:length(real_eig1p)
        for NN=1:nDim
            G_RP(M,N)=G_RP(M,N)+ (a_mat(M,NN)*A(NN) - AUX_mat(M,NN)*sqrt(B1(N+min(real_Ray1_idx)-m1)))*FRP(NN,N);
        end
    end
end

for M=1:nDim
    for N=1:length(real_eig2n)
        for NN=1:nDim
           G_RN(M,N)=G_RN(M,N)+ (a_mat(M,NN)*A(NN) + AUX_mat(M,NN)*sqrt(B2(N+min(real_Ray2_idx)-m1)))*FRN(NN,N);
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
%clearvars -except etaR etaT Lam VecWaveLength Parameters d h thI n1 n2 K mu0 eps0 nTr tol alfa0 mu eps1 eps2 mu1 mu2 nDim real_Ray1_idx i RefVec
etaT;
etaR;
sum(etaT);
sum(etaR);
RAmp=(real_Ray1_idx==0);
RefVec(i)=etaR(RAmp);
clear real_Ray1_idx
toc
end
plot(VecWaveLength,RefVec);
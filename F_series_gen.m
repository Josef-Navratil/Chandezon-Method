%% Generates Fouries series for given periodic function
function[f_vec] =  F_series_gen(a_fun,m,T,nDim) 
%input parameters - function profile, order of fft, period
%length, number of terms in the series
N=2^m;
tol=10e-10;
cut_small=1;
y_vec=linspace(0,T,N+1); 
f = a_fun(y_vec);
fhat =f ;
fhat(1)=(f(1)+f(N+1))/2;
fhat(N+1)=[];

F=fft(fhat,N)/N;
if (cut_small==1)
    ind_small_real=(abs(real(F))<tol);
    F(ind_small_real)=1i.*imag(F(ind_small_real));
    ind_small_imag=(abs(imag(F))<tol);
    F(ind_small_imag)=real(F(ind_small_imag));
    f_vec=F(1:nDim)';
end
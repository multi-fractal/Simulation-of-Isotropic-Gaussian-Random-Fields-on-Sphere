
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models parameters

NE=100000;   % Number of realizations

Nh=100;      % Number of harmonics

% First point on sphere
phi1=pi/6;
theta1=pi/6;

% Second point on sphere
% phi2=pi/3;
% theta2=pi/3;
phi2=3*pi/2;
theta2=pi/2;
% phi2=pi/6;
% theta2=pi/6;

x1=cos(phi1)*sin(theta1);
y1=sin(phi1)*sin(theta1);
z1=cos(theta1);

x2=cos(phi2)*sin(theta2);
y2=sin(phi2)*sin(theta2);
z2=cos(theta2);

alpha=acos(x1*x2+y1*y2+z1*z2);

Nmax=51; % spectrum truncation

a=zeros(1,Nmax);

% Spectrum

z=0.8;

for iN=1:Nmax
    
  a(iN) = z^(iN-1); 
    
end;

pa=a/sum(a);

suma=sum(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of spherical harmonics

Ynk1=zeros(Nmax,2*Nmax-1);
Ynk2=zeros(Nmax,2*Nmax-1);

for iN=1:Nmax
    
N=iN-1;

    Pnk10=legendre(N,cos(theta1));
    Pnk20=legendre(N,cos(theta2));
    
    for i=1:N+1
        
        Cnk= sqrt( (0.25*(2*N+1)/pi) * factorial(N-(i-1))/factorial(N+(i-1)) );        
        Ynk1(iN,i+N)=Cnk*Pnk10(i)*exp(1i*(i-1)*phi1);
        Ynk2(iN,i+N)=Cnk*Pnk20(i)*exp(1i*(i-1)*phi2);
        
    end;
    
    for i=1:N
        
        Ynk1(iN,N-i+1)=(-1)^(i)*conj(Ynk1(iN,i+N+1));
        Ynk2(iN,N-i+1)=(-1)^(i)*conj(Ynk2(iN,i+N+1));
        
    end; 

end;    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ModeliNh

C=0;

for iNE=1:NE

iNE    

T1=0;
T2=0;

for ig=1:Nh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    sumP=0;
    r=rand;
    N=-1;
    
    while (sumP<r)
        
        N=N+1;
        sumP=sumP+pa(N+1);
        
    end;
    
    iN=N+1;

    K=floor(rand*(2*N+1))+1;
    
    ksi=randn;
    eta=randn;

    T1=T1+2*sqrt(suma*pi) * ( ksi*real(Ynk1(iN,K)) + eta*imag(Ynk1(iN,K)) );
    T2=T2+2*sqrt(suma*pi) * ( ksi*real(Ynk2(iN,K)) + eta*imag(Ynk2(iN,K)) );   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

end;

T1=T1/sqrt(Nh);
T2=T2/sqrt(Nh);

C=C+T1*T2;

end;

% Analytical covariance
C_an=1/sqrt(1-2*z*cos(alpha)+z^2)
% Numerically estimated covariance
C=C/NE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

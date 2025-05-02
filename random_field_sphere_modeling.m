
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates for drawing on a sphere

Nx=128;
Ny=Nx;

[X, Y, Z]=sphere(Nx-1);
[AZ, EL, R] = cart2sph(X, Y, Z);

AZ=AZ+pi;
EL=EL+pi/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrum

Nmax=51; % spectrum truncation

a=zeros(1,Nmax);

Za=0.8;

for iN=1:Nmax
    
  a(iN) = Za^(iN-1); 
    
end;

pa=a/sum(a);

suma=sum(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of spherical harmonics

Ynk=zeros(Nmax,2*Nmax-1,Nx,Ny);

for ix=1:Nx
for iy=1:Ny

    phi=AZ(ix,iy);    
    theta=EL(ix,iy);    
    
    for iN=1:Nmax
    
    N=iN-1;

        Pnk0=legendre(N,cos(theta));
    
        for i=1:N+1
        
            Cnk= sqrt( (0.25*(2*N+1)/pi) * factorial(N-(i-1))/factorial(N+(i-1)) );        
            Ynk(iN,i+N,ix,iy)=Cnk*Pnk0(i)*exp(1i*(i-1)*phi);
        
        end;
    
        for i=1:N
        
            Ynk(iN,N-i+1,ix,iy)=(-1)^(i)*conj(Ynk(iN,i+N+1,ix,iy));
        
        end; 

    end;   

end;
end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ModeliNh

Nh=1000;         % Number of harmonics

T=zeros(Nx,Ny);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for ig=1:Nh        
      
        ksi=randn; 
        eta=randn; 
        
        Tpt=zeros(Nx,Ny);
    
        sumP=0;
        r=rand;
        N=-1;
    
        while (sumP<r)
        
            N=N+1;
            sumP=sumP+pa(N+1);
        
        end;
    
        iN=N+1;

        K=floor(rand*(2*N+1))+1;

        for ix=1:Nx
        for iy=1:Ny

            Tpt(ix,iy)=Tpt(ix,iy)+2*sqrt(suma*pi) * ( ksi*real(Ynk(iN,K,ix,iy)) + eta*imag(Ynk(iN,K,ix,iy)) );
        
        end;
        end;
           
        T=T+Tpt;
        
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=T/sqrt(Nh);   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
figure;
surf(X,Y,Z,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
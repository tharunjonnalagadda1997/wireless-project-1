clc;
clear all;
close all;
s=0;
Kc=10;  % no of CUs
Kp=10;  % no of PUs 
N=100;  % no of CBS antennas
N_iter=100; % no of iterations
N_iter1=100;
sp=0.8; % correlation coefficient
Pr=[];
geC=zeros(N,1);
geC=zeros(N,1);
heC=zeros(N,1);
yp=zeros(N_iter,1);
T=zeros(N_iter,1);
B=zeros(N_iter,1);
Q1=zeros(N_iter1,1);
for eta=1:0.1:3   
for p1=1:1:N_iter1
Alpha=10^(-6)+(10^(-3)-10^(-6))*rand(Kp,1);
Gamma=10^(-6)+(10^(-3)-10^(-6))*rand(Kc,1);
Phi=10^(-6)+(10^(-3)-10^(-6))*rand(Kc,1);  
for p=1:1:N_iter
% declaration of estimate channel vector between PUs and CBS antennas
for j=1:1:Kp
geC(:,j)=sqrt(Alpha(j)/2).*(randn(N,1)+1i*randn(N,1));
end
% declaration of etimate error channel vector between PUs and CBS antennas
for j=1:1:Kp
egC(:,j)=sqrt(Alpha(j)/2).*(randn(N,1)+1i*randn(N,1));
end
% declaration of true channel vector between PUs and CBS antennas
gC=sp*geC(:,:)+sqrt(1-sp^2)*egC(:,:);
% declaration of estimate channel vector between CUs and CBS antennas
for j=1:1:Kc
heC(:,j)=sqrt(Gamma(j)/2).*(randn(N,1)+1i*randn(N,1));
end
% declaration of true powers for normalised precoding vector
for i=1:1:Kp
U=0;
for j=1:1:Kc
U=U+(1/Kc).*abs(transpose(gC(:,i))*conj(heC(:,j)/norm(heC(:,j))))^2;
end
Z(i,p)=U;
end
% declaration of estimate powers for normalised precoding vector
for i=1:1:Kp
W=0;
for j=1:1:Kc
W=W+(1/Kc).*abs(transpose(geC(:,i))*conj(heC(:,j)/norm(heC(:,j))))^2;
end
Ze(i,p)=W;
end
% finding out max true power among all PUs for given iteration
T(p)=max(Z(:,p));
% finding out max estimate power among all PUs for given iteration
B(p)=max(Ze(:,p));
% checking the inequality condition for exact equality case for normalised
% precoding vector
if (T(p)<=eta*B(p))
l=1;
else
l=0;
end
yp(p)=l;
end
Q=sum(yp(:))/N_iter;
Q1(p1)=Q;
end
Q3=sum(Q1(:,1))/N_iter1;
Pr=[Pr Q3];
end
eta=1:0.1:3;
%figure;
plot(eta,Pr,'*--');
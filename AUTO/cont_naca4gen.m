% Continuous NACA 4 gen

function [coordU,coordL] = cont_naca4gen(x,m,p,t)

a0= 0.2969;
a1=-0.1260;
a2=-0.3516;
a3= 0.2843;
a4=-0.1036;  % For zero thick TE

x = x';
yt=(t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

I1 = x <= p;
I2 = x > p;
xc1=x(find(x<=p));
xc2=x(find(x>p));
xc=[xc1 ; xc2];

if p==0
    xu=x;
    yu=yt;

    xl=x;
    yl=-yt;
    
    zc=zeros(size(xc));
else
    yc1=(m/p^2)*(2*p*xc1-xc1.^2);
    yc2=(m/(1-p)^2)*((1-2*p)+2*p*xc2-xc2.^2);
    zc=[yc1 ; yc2];

    dyc1_dx=(m/p^2)*(2*p-2*xc1);
    dyc2_dx=(m/(1-p)^2)*(2*p-2*xc2);
    dyc_dx=[dyc1_dx ; dyc2_dx];
    theta=atan(dyc_dx);

    xu=x-yt.*sin(theta);
    yu=zc+yt.*cos(theta);

    xl=x+yt.*sin(theta);
    yl=zc-yt.*cos(theta);
end

coordU = [xu, yu]';
coordL = [xl, yl]';
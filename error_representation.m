syms x y z v 
format long g

%%%%%% Analysis of accuracy of computer computation

% main function z
z = ((x) + (cos(y)/3))/((y^2) - (sin(y)/2));

% using functions of error propagation
% that are used to check if my 
% calculations by hand are correct
Tx=x/z*diff(z,x);
Ty=y/z*diff(z,y);

zsin=subs(z,sin(y),v);
Ksin=v/zsin*diff(zsin,v);
Ksin=subs(Ksin,v,sin(y));

zcos=subs(z,cos(y),v);
Kcos=v/zcos*diff(zcos,v);
Kcos=subs(Kcos,v,cos(y));

zby2=subs(z,sin(y)/2,v);
Kby2=v/zby2*diff(zby2,v);
Kby2=subs(Kby2,v,sin(y)/2);

zby3=subs(z,cos(y)/3,v);
Kby3=v/zby3*diff(zby3,v);
Kby3=subs(Kby3,v,cos(y)/3);

zpow=subs(z,y^2,v);
Kpow=v/zpow*diff(zpow,v);
Kpow=subs(Kpow,v,y^2);

zadd=subs(z,x+cos(y)/3,v);
Kadd=v/zadd*diff(zadd,v);
Kadd=subs(Kadd,v,x+cos(y)/3);

zsub=subs(z,y^2-sin(y)/2,v);
Ksub=v/zsub*diff(zsub,v);
Ksub=subs(Ksub,v,y^2-sin(y)/2);

zdiv=subs(z,((x) + (cos(y)/3))/((y^2) - (sin(y)/2)),v);
Kdiv=v/zdiv*diff(zdiv,v);
Kdiv=subs(Kdiv,v,((x) + (cos(y)/3))/((y^2) - (sin(y)/2)));

matlab_diff = [ Tx Ty Ksin Kcos Kby2 Kby3 Kpow Kadd Ksub Kdiv ];

%% Task 1
% formulas obtained from analytical differentiation by hand
TAx=(x/(x+(cos(y)/3)));
TAy=(y/(x+(cos(y)/3)))* ...
    (((-(sin(y))/3)*((y^2)-((sin(y))/2))-((((cos(y))/3)+x) ...
    *(2*y-((cos(y))/2))))/((y^2)-((sin(y))/2)));
KAsin=(sin(y)/(2*y^2-sin(y)));
KAcos=(cos(y)/(3*x+(cos(y)/3)));
KAby2=(sin(y)/(2*y^2-sin(y)));
KAby3=(cos(y)/(3*x+(cos(y)/3)));
KApow=(-(y^2)/(y^2-(sin(y)/2)));
KAdiv=1;
KAsub=-1;
KAadd=1;

hand_diff = [ TAx TAy KAsin KAcos KAby2 KAby3 KApow KAadd KAsub KAdiv ];

% formulas obtained from epsilon calculus by hand
TEx=(x/(x+(cos(y)/3)));
TEy=-((2*y*sin(y)*((y^2)-(sin(y)/2)))+3*(4*(y^2)-y*cos(y)) ...
    *(((cos(y))/3)+x))/(6*((((cos(y))/3)+x)*(((y^2)-(sin(y)/2)))));
KEsin=(sin(y)/(2*y^2-sin(y)));
KEcos=(cos(y)/(3*x+(cos(y)/3)));
KEby2=(sin(y)/(2*y^2-sin(y)));
KEby3=(cos(y)/(3*x+(cos(y)/3)));
KEpow=-((y^2)/(y^2-(sin(y)/2)));
KEdiv=1;
KEsub=-1;
KEadd=1;

eps_calc = [ TEx TEy KEsin KEcos KEby2 KEby3 KEpow KEadd KEsub KEdiv ];

names = [ "Tx" "Ty" "Ksin" "Kcos" "Kby2" ...
                "Kby3" "Kpow" "Kadd" "Ksub" "Kdiv" ];

for n = 1 : length(eps_calc)
    figure(n+10);
    subplot(1,2,1)
    fsurf(eps_calc(n) , [0,1,0,1])
    title("Analytical differentiation by hand");
    xlabel('x');
    ylabel('y');
    subplot(1,2,2)
    fsurf(hand_diff(n) , [0,1,0,1])
    title("Epsilon calculus");
    xlabel('x');
    ylabel('y');
    sgtitle(names(n));
end

%% Task 2
%%%%% indicator of the accuracy of the floating-point representation is eps ,
%%%%% assess the total error of computing the value of z by maximising 
%%%%% indicator d1
eps = 2*10^(-12);

d1=abs(TAx)+abs(TAy)+abs(KAsin)+abs(KAcos)+abs(KAby2)+abs(KAby3) ...
    +abs(KApow)+abs(KAdiv)+abs(KAsub)+abs(KAadd);

d2 = matlabFunction(d1);
fsurf(d2, [1,10,1,10]);
d3 = max(d2(x1,y1))*eps;

maximal = double(d3);

%% Task 3 
%%%%% Compare the result obtained in task 2 with the estimation
%%%%% obtained by means of the simulation method
% initialize matrices 
matrix1=de2bi(0:1023);
matrix2=ones(size(matrix1))*eps;

sol=sym(ones(1,1024));
x1 = linspace(1,10,100);
y1 = linspace(1,10,100);

for i=1:1024
    z0=TAx*matrix2(i,1)+TAy*matrix2(i,2)+KAsin*matrix2(i,3) ...
    +KApow*matrix2(i,4)+KAcos*matrix2(i,5)+KAadd*matrix2(i,6) ...
    +KAsub*matrix2(i,7)+KAdiv*matrix2(i,8)+KAby2*matrix2(i,9) ...
    +KAby3*matrix2(i,10);
    z4 = matlabFunction(abs(z0));
    sol(1,i) = max(z4(x1,y1));
end 

maximum = double(max(sol));

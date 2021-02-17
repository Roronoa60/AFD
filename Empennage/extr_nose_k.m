%This is a function for obtaining the K values for curved plates in shear. Inputs are a/b and b/sqrt(Rt) as specified on
%the ESDU 02.03.18 graph.

%Note that a=a/b and b= b/sqrt(Rt)

function K=extr_nose_k(b,a)

load('d_sect.mat');

if b>20 || b<0
    disp('Value of b/sqrt(Rt) is out of accurate range. Beware!')
elseif a<1 || (a>3 && a~=inf)
    disp('Value of a/b is out of accurate range. Beware!')
end 
    
switch a
    case 1
        K=feval(FIT_1,b);
    case 1.5
        K=feval(FIT_2,b);
    case 2
        K=feval(FIT_3,b);
    case 3
        K=feval(FIT_4,b);
    case inf
        K=feval(FIT_5,b);
    otherwise %Surface fit. Note that a/b=10 is used as the inf case for the surface fit.
        K=FIT_TOTAL(b,a)
end
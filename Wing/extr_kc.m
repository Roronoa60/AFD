%This is a function for obtaining the K_s values for shear for buckling
%of flat plate under in-plane bending loads. Input is a/b as specified on
%the Niu graph. Specify type 1,2,3,4 according to boundary conditions.

function K=extr_kc(a,type)

load('kc_var.mat');

if a>6 || a<0
    disp('Value of a/b is out of accurate range. Beware!')
end 
    
switch type
    case 1
        K=feval(FIT_1,a);
    case 2
        K=feval(FIT_2,a);
    case 3
        K=feval(FIT_3,a);
    case 4
        K=feval(FIT_4,a);
    case 5
        K=feval(FIT_5,a);
    case 6
        K=feval(FIT_6,a);
    case 7
        K=feval(FIT_7,a);
    case 8
        K=feval(FIT_8,a);
    case 9
        K=feval(FIT_9,a);
    case 10
        K=feval(FIT_10,a);
    otherwise
        disp('Please Enter Correct Constraint Mode')
end
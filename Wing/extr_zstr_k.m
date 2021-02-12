%This is a function for obtaining the K values for z-stringers on a stringer-skin panel. Input is h/b as specified on
%the ESDU graphs. Specify d/h=0.3,0.4,0.5.

function K=extr_zstr_k(a,t_rat,d_h)

load('kz_dh3.mat');
load('kz_dh4.mat');
load('kz_dh5.mat');

%t_rat = ts/t
%d_h=d/h=0.3 or 0.4 or 0.5
d_h_checker=[0.3 0.4 0.5];

if a>1 || a<0.1
    disp('Value of h/b is out of accurate range. Beware!')
end 
    
if t_rat<0.5 || t_rat>2
    disp('Value of ts/t is out of accurate range. Beware!')
end 

% Interpolation for if t_rat is not discrete number
if (0.5 <= t_rat) & (t_rat <= 1)
   t_rat = round(t_rat, 1);
elseif t_rat <= 1.125
    t_rat = 1;
elseif t_rat <= 1.375
    t_rat = 1.25;
elseif t_rat <= 1.75
    t_rat = 1.5;
elseif t_rat <= 3
    t_rat = 2;
end

if any(d_h_checker(:) == d_h)~=1
    error('d/h values other than 0.3, 0.4, 0.5 are not supported')
end 

%% Create massive conditional statement for d/h values

if d_h==0.3
    %Now if t_rat actually is one of the curves provided
    switch t_rat
        case 0.5
            K=FIT_1_3(a);
        case 0.6
            K=FIT_2_3(a);
        case 0.7
            K=FIT_3_3(a);
        case 0.8
            K=FIT_4_3(a);
        case 0.9
            K=FIT_5_3(a);
        case 1
            K=FIT_6_3(a);
        case 1.25
            K=FIT_7_3(a);
        case 1.5
            K=FIT_8_3(a);
        case 2
            K=FIT_9_3(a);
        otherwise
            K=FIT_TOTAL_3(a,t_rat);       
    end 
        
elseif d_h==0.4
    
        %Now if t_rat actually is one of the curves provided
    switch t_rat
        case 0.5
            K=FIT_1_4(a);
        case 0.6
            K=FIT_2_4(a);
        case 0.7
            K=FIT_3_4(a);
        case 0.8
            K=FIT_4_4(a);
        case 0.9
            K=FIT_5_4(a);
        case 1
            K=FIT_6_4(a);
        case 1.25
            K=FIT_7_4(a);
        case 1.5
            K=FIT_8_4(a);
        case 2
            K=FIT_9_4(a);
        otherwise
            K=FIT_TOTAL_4(a,t_rat);       
    end  
    
else %d_h==0.5
        %Now if t_rat actually is one of the curves provided
    switch t_rat
        case 0.5
            K=FIT_1_5(a);
        case 0.6
            K=FIT_2_5(a);
        case 0.7
            K=FIT_3_5(a);
        case 0.8
            K=FIT_4_5(a);
        case 0.9
            K=FIT_5_5(a);
        case 1
            K=FIT_6_5(a);
        case 1.25
            K=FIT_7_5(a);
        case 1.5
            K=FIT_8_5(a);
        case 2
            K=FIT_9_5(a);
        otherwise
            K=FIT_TOTAL_5(a,t_rat);       
    end 
end 

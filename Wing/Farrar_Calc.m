function F = Farrar_Calc(a, t_rat)
%This fucntion calculates the Farrar efficiency factor for a given ts/t =
%t_rat and As/bt = a.

load('FARRAR.mat')

K=FIT(a,t_rat)

end 
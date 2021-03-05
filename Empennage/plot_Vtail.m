close all

t_blue = [0.3 0.6 1];
orange = [1, 0.6, 0];
green = [0.3 0.7 0.2];
red = [0.9 0.1 0.2];

L = extract_dimension(0, 'wing_span', 'VT');
C0 = extract_dimension(0, 'chord', 'VT');
Ct = extract_dimension(L, 'chord', 'VT');
sweep = 29.8; % HT Same
b1 = 0.143;
rib_position = [0,0.503885447033903,1.01793040648245,1.54114716947606,2.07441751493461,2.61934611178983,3.17720823904346,3.74588733981849,4.32197785484113,4.90494310410758,5.49880429310460,6.20336711995726,7];
draw_wing(C0 ,Ct ,L ,sweep ,b1 , rib_position)

function draw_wing(C0 ,Ct ,L ,sweep ,b1 , rib_position)
% draw_wing(C0,Ct,L,sweep, b1, rib_position) plot wing geometry with
% constant stringer pitch, b1 and rib spacing
% INPUT:
% C0 - chord length at root [m]
% Ct - chord length at tip [m]
% L - half span [m]
% sweep - sweep angle of Leading edge [m]
% b1 - pitch of stringer [m]
% rib_position - positon of the rib (array) [m]
% OUTPUT:
% figure plot of the wing
fs = 0.18 ; % change to 0.15 for htail  0.18 for vtail
bs = 0.7 ; % change to 0.7 for htail  0.7 for vtail
tsa = tand(sweep) ; % tan of sweep angle of leading egde
fsa = atand( (L*tsa-C0*fs+Ct*fs) / L ); %angle of front spar
bsa = atand( (L*tsa-C0*bs+Ct*bs) / L ) ;
ssa= fsa; % stringer sweep angle = front spar sweep angle

figure
set(gca,'TickDir','out'); 
ylabel('Spanwise Location (m)')
xlabel('Chordwise Location (m)')
axis equal
% plot wing with both front and rear spar
line([0 C0 (L*tsa+Ct) (L*tsa) 0],[0 0 L L 0], 'LineWidth',2,'Color','k');hold on
line([(C0*bs) (L*tsa+Ct*bs) (L*tsa+Ct*fs) (C0*fs)], [0 L L 0],'LineWidth',2,'Color','b')

% plot stringer
if ssa >= bsa
    for d = C0*fs +b1 : b1 : C0*bs % back spar
        x_sp = (C0*bs - d) / (tand(ssa)-tand(bsa)) ;
        y_sp = tand(ssa)*(x_sp) + d ;
        if y_sp <= tand(fsa)*(x_sp) + C0*fs
            y_sp = tand(fsa)*(x_sp) + C0*fs;
        end
        line([d y_sp],[0 x_sp],'Color','r', 'Linewidth', 0.05)
    end
else
    for d = C0*fs : b1 : C0*bs % front spar
        x_sp = (C0*fs - d) / (tand(ssa)-tand(fsa)) ;
        y_sp = tand(ssa)*(x_sp) + d ;
        if y_sp >= tand(bsa)*(x_sp) + C0*bs
            y_sp = tand(bsa)*(x_sp) + C0*bs;
        end
        line([d y_sp],[0 x_sp],'Color','r')
    end
end

% plot the ribs
c1 = tand(fsa).* rib_position + C0*fs ;
x_rib= (-C0*bs- tand(90+bsa).*rib_position + c1) ./ (tand(bsa) - tand(90+bsa)) ;
y_rib= tand(90+bsa).*(x_rib - rib_position) + c1 ;
for i = 1:size(c1,2)
    line([y_rib(i) c1(i)],[x_rib(i) rib_position(i)],'LineWidth',1.8,'Color','b')
end

% fixed the axis
axis([0 L*tsa+Ct 0 L])
end
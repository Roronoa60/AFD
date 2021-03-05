close all
clear
clc
addpath('Stringer Panel Result Data')

t_blue = [0.3 0.6 1];
orange = [1, 0.6, 0];
green = [0.3 0.7 0.2];
red = [0.9 0.1 0.2];
select_color = [0 0 0];
select_size = 60;
select_thickness = 2.2;

load('bottom_weights')
bottom_n = cell2mat(bottom_weights(:, 1));
bottom_weight = cell2mat(bottom_weights(:, 4));

line([30 30, 0],[0 2.7*bottom_weight(13), 2.7*bottom_weight(13)], 'LineStyle','--', 'Color', 'black', 'LineWidth',1.5)
hold on
%scatter(bg.n, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',1.8)
%scatter(bottom_n, 2.7*bottom_weight(13),'x', 'MarkerEdgeColor', red, 'LineWidth',0.2)
%scatter(30, 2.7*bottom_weight(13), select_size, 'MarkerEdgeColor', select_color, 'LineWidth',1.8)
scatter(bottom_n, 2.7*bottom_weight, 80,'x', 'MarkerEdgeColor', green, 'LineWidth',2)
line([30 30, 0],[0 2.7*bottom_weight(13), 2.7*bottom_weight(13)], 'LineStyle','--', 'Color', 'black', 'LineWidth',1.5)
xlabel('Number of Stringers - n');
ylabel('Weight (kg)');
legend('Selected Value', 'Best Weight Configuration')
ylim([0 1600])

figure

addpath('Stringer Panel Result Data')
load('weights_n50.mat')
%load('weights_n50_best_geometry.mat')
load('n30_best_geometry')
bg = best_geometry;
titles = {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'L', 'root weight', 'total weight', 'As_bt', 'te'};

weights = weights_array;
n = cell2mat(weights(:, 1));
t = cell2mat(weights(:, 2));
Rt = cell2mat(weights(:, 3));
Rdh = cell2mat(weights(:, 4));
Rb = cell2mat(weights(:, 5));
F = cell2mat(weights(:, 6));
L = cell2mat(weights(:, 7));
root_weight = cell2mat(weights(:, 8));
total_weight = cell2mat(weights(:, 9))*2.7;
As_bt = cell2mat(weights(:, 10));
te = cell2mat(weights(:,14));
T = table(n, t, Rt, Rdh, Rb, F, L, root_weight, total_weight, As_bt, te);


subplot(2, 2, 1)
scatter(bg.F, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',select_thickness)
hold on

for i=8:2:50
    current_n = T(T.n==i, :);
    scatter(current_n.F, current_n.total_weight,'.', 'MarkerEdgeColor', t_blue);
    set(gca,'FontSize',12)
end
scatter(bg.F, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',select_thickness)
legend('Selected Value')
xlabel('FARRAR Efficiency Factor');
ylabel('Weight (kg)');
ylim([0 2500])

subplot(2, 2, 2)
scatter(bg.L_distribution(2), bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',select_thickness)
hold on
for i=8:2:50
    current_n = T(T.n==i, :);
    scatter(current_n.L, current_n.total_weight, '.', 'MarkerEdgeColor', t_blue');
    
    set(gca,'FontSize',12)
end
scatter(bg.L_distribution(2), bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',select_thickness)
legend('Selected Value')
xlabel('Average Rib Spacing (m)');
ylabel('Weight (kg)');
ylim([0 2500])
xlim([0 1])

subplot(2, 2, 3)
scatter(bg.As_bt, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',2)
hold on
for i=8:2:50
    current_n = T(T.n==i, :);
    scatter(current_n.As_bt, current_n.total_weight,'.', 'MarkerEdgeColor', t_blue);
    set(gca,'FontSize',12)
end
scatter(bg.As_bt, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',select_thickness)
legend('Selected Value')
xlabel('As/bt');
ylabel('Weight (kg)');
ylim([0 2500])

sorted_weights = sortrows(T, 'total_weight', 'ascend');



Ribs_weight = 279.72;
subplot(2, 2, 4)
hold on
line([30 30, 0],[0 bg.total_weight*2.7, bg.total_weight*2.7], 'LineStyle','--', 'Color', 'black', 'LineWidth',1.5)
%scatter(bg.n, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',1.8)
scatter(n, total_weight,'x', 'MarkerEdgeColor', red, 'LineWidth',0.2)
for i=8:2:50
    current_n = T(T.n==i, :);
    max_n = T(T.total_weight == min(current_n.total_weight), :);

    weight_ratio = max_n.As_bt;
    
    scatter(max_n.n, max_n.total_weight, 80,'x', 'MarkerEdgeColor', green, 'LineWidth',2);
    %scatter(max_n.n, (max_n.total_weight-Ribs_weight)*weight_ratio*(1+0.002*i),'x', 'MarkerEdgeColor', orange, 'LineWidth',1.5);
    %scatter(max_n.n, (max_n.total_weight-Ribs_weight)*(1-weight_ratio*(1+0.002*i)),'x', 'MarkerEdgeColor', green, 'LineWidth',1.5);
    set(gca,'FontSize',12)
end
line([30 30, 0],[0 bg.total_weight*2.7, bg.total_weight*2.7], 'LineStyle','--', 'Color', 'black', 'LineWidth',1.5)
%scatter(bg.n, bg.total_weight*2.7, select_size, 'MarkerEdgeColor', select_color, 'LineWidth',1.8)
legend('Selected Value', 'Suboptimal Configuration', 'Best Weight Configuration')
%line([30 30, 0],[0 1248.82, 1248.82], 'LineStyle','--', 'Color', 'black', 'LineWidth',1.5)

ylim([0 2500])
xlim([0 50])

ylabel('Weight (kg)');
xlabel('Number of Stringers - n');

figure
scatter(bg.L_distribution(2:end-1), bg.Tr_distribution(1:end-1), 'LineWidth',1.5, 'MarkerEdgeColor', t_blue);
hold on
stairs(bg.L_distribution, [bg.t_distribution(1:end-1), bg.t_distribution(end-1)], 'LineWidth',1.5, 'Color', red);
stairs(bg.L_distribution(1:end), bg.bottom_t, 'LineWidth',1.5, 'Color', green);

xlabel('Wing Span (m)')
ylabel('Thickness (mm)')
ylim([0 10]);
xlim([0 18]);
set(gca,'FontSize',12)
legend('Rib Thickness', 'Top Skin Thickness', 'Bottom Skin Thickness')

total_weight = 2700*(bg.total_weight + bg.bottom_weight)/1000

stringer_weight_ratio = max_n.As_bt * 30/(30+1);

top_stringer_weight = stringer_weight_ratio*969.03
bottom_stringer_weight = bg.bottom_weight*2.7*stringer_weight_ratio
top_skin_weight = (1-stringer_weight_ratio)*969.03

bottom_skin = bg.bottom_weight*2.7 - top_stringer_weight

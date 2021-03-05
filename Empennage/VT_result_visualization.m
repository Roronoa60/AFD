close all
clear
clc

t_blue = [0.3 0.6 1];
red = [0.9 0.1 0.2];
green = [0.3 0.7 0.2];

addpath('Results Data')
%load('varying_L_weights_array.mat')
%load('varying_L_best_geometry.mat')
%load('bottom_stringer_for_n40')
load('VT_weights_n40_V4.mat')
load('best_geometry_V4.mat')

titles = {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'L', 'root weight', 'total weight', 'As_bt', 'te'};

%weights = cell2mat(weights);
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

%{
subplot(2, 2, 1)
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.F, current_n.total_weight,'.', 'MarkerEdgeColor', t_blue);
    hold on
end
%legend(string(4:2:40));
xlabel('FARRAR Efficiency Factor');
ylabel('Weight (kg)');
ylim([0 450])

%{
figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.root_weight, current_n.total_weight);
    hold on
end
legend(string(4:2:40));
xlabel('root weight');
ylabel('total weight');
%}
%scatter(root_weight, total_weight);

subplot(2, 2, 2)
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.L, current_n.total_weight, '.', 'MarkerEdgeColor', t_blue');
    hold on
end
%legend(string(4:2:40));
xlabel('Average Rib Spacing (m)');
ylabel('Weight (kg)');
ylim([0 450])
xlim([0 3])



%{
figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.te, current_n.total_weight);
    hold on
end

legend(string(4:2:40));
xlabel('te');
ylabel('total weight');
%}

subplot(2, 2, 3)
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.As_bt, current_n.total_weight,'.', 'MarkerEdgeColor', t_blue);
    hold on
end

%legend(string([4, 6, 10, 16, 24, 32, 40]));
xlabel('As/bt');
ylabel('Weight (kg)');
ylim([0 450])

sorted_weights = sortrows(T, 'total_weight', 'ascend');
%}
%subplot(2, 2, 4)
figure
set(gca,'FontSize',12)
line([20 20, 0],[0 185.6, 185.6], 'LineStyle','--', 'Color', 'black', 'LineWidth',1.5)
hold on
scatter(n, total_weight,'x', 'MarkerEdgeColor', red, 'LineWidth',0.2)

for i=4:2:40
    current_n = T(T.n==i, :);
    max_n = T(T.total_weight == min(current_n.total_weight), :);

    weight_ratio = max_n.As_bt;
    
    scatter(max_n.n, max_n.total_weight, 80,'x', 'MarkerEdgeColor', green, 'LineWidth',2);
    %scatter(max_n.n, max_n.total_weight*weight_ratio,'x', 'MarkerEdgeColor', red, 'LineWidth',1.5);
    %scatter(max_n.n, max_n.total_weight*(1-weight_ratio),'x', 'MarkerEdgeColor', green, 'LineWidth',1.5);
    
end
legend('Selected Value', 'Suboptimal Configuration', 'Best Weight Configuration') %, 'Combined Weight', 'Stringer Weight Component', 'Skin Weight Component')

xlabel('Number of Stringers - n');
ylim([0 450])
ylabel('Weight (kg)');




load('VT_n20_best_geometry')
load('VT_D_section')
D_t_distribution = cell2mat(D_section(4));

bg = best_geometry;
figure
scatter(bg.L_distribution(2:end), bg.Tr_distribution(1:end), 'LineWidth',1.5, 'MarkerEdgeColor', red);
set(gca,'FontSize',12)
hold on
stairs([bg.L_distribution, extract_dimension(0, 'wing_span', 'VT')], [bg.t_distribution, bg.t_distribution(end)], 'LineWidth',1.5, 'Color', green);
stairs([bg.L_distribution, extract_dimension(0, 'wing_span', 'VT')], [D_t_distribution, D_t_distribution(end)], 'LineWidth',1.5, 'Color', t_blue);
xlabel('Wing Span (m)')
ylabel('Thickness (mm)')
ylim([0 10]);

legend('Rib and Pseudo Rib Thickness', 'Skin Thickness', 'D Section Thickness')

%scatter(bg.L_distribution(1:end), best_n.t_distribution, 'x');
%legend('Rib Thickness', 'Top Skin Thickness', 'Bottom Skin Thickness')
%}
%{
current_n = T(T.20==i, :);
chosen = T(T.total_weight == min(current_n.total_weight), :)
scatter(chosen.L_distribution(2:end), bg.Tr_distribution)
%}
total_weight = 2700*(bg.total_weight/1000 + cell2mat(D_section(3)))
D_section_weight = cell2mat(D_section(3))*2700
rib_weight = 0.138*(total_weight-D_section_weight)
stringer_weight = (bg.As_bt*20/21)*(total_weight-D_section_weight-rib_weight)
skin_weight = total_weight-D_section_weight-rib_weight-stringer_weight
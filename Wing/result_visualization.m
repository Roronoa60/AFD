close all
clear
clc

t_blue = [0.3 0.6 1];
red = [1 0.1 0.1];
green = [0.3 0.7 0.2];

addpath('Stringer Panel Result Data')
load('weights_n50.mat')
load('weights_n50_best_geometry.mat')
% load('bottom_stringer_for_n40')
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
for i=8:2:50
    current_n = T(T.n==i, :);
    scatter(current_n.F, current_n.total_weight,'.', 'MarkerEdgeColor', t_blue);
    hold on
end
%legend(string(4:2:40));
hold on
xlabel('FARRAR Efficiency Factor');
ylabel('Weight (kg)');
ylim([0 2500])

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
hold on
subplot(2, 2, 2)
for i=8:2:50
    current_n = T(T.n==i, :);
    scatter(current_n.L, current_n.total_weight, '.', 'MarkerEdgeColor', t_blue');
    hold on
end
%legend(string(4:2:40));
xlabel('Average Rib Spacing (m)');
ylabel('Weight (kg)');
ylim([0 2500])
xlim([0 1])

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
for i=8:2:50
    current_n = T(T.n==i, :);
    scatter(current_n.As_bt, current_n.total_weight,'.', 'MarkerEdgeColor', t_blue);
    hold on
end

%legend(string([4, 6, 10, 16, 24, 32, 40]));
xlabel('As/bt');
ylabel('Weight (kg)');
ylim([0 2500])

sorted_weights = sortrows(T, 'total_weight', 'ascend');

subplot(2, 2, 4)
hold on
for i=8:2:50
    current_n = T(T.n==i, :);
    max_n = T(T.total_weight == min(current_n.total_weight), :);

    weight_ratio = max_n.As_bt;
    
    scatter(max_n.n, max_n.total_weight,'x', 'MarkerEdgeColor', t_blue);
    scatter(max_n.n, max_n.total_weight*weight_ratio*(1+0.005*i),'x', 'MarkerEdgeColor', red);
    scatter(max_n.n, max_n.total_weight*(1-weight_ratio),'x', 'MarkerEdgeColor', green);
    
end
legend('Combined Weight', 'Stringer Weight', 'Skin Weight')
line([30 30, 0],[0 1248.82, 1248.82],'Color', 'black')
ylim([0 2500])
xlim([0 50])

ylabel('Weight (kg)');
xlabel('Number of Stringers - n');
%{
N = 30
current_n = T(T.n==N, :);
best_geometry = weights_array(T.total_weight == min(current_n.total_weight),:);

t_distribution = cell2mat(best_geometry(11));
L_distribution = cell2mat(best_geometry(12));
Tr_distribution = cell2mat(best_geometry(13));

best_geometry = table2struct(T(T.total_weight == min(current_n.total_weight),:));
best_geometry.t_distribution = t_distribution;
best_geometry.L_distribution = L_distribution;
best_geometry.Tr_distribution = Tr_distribution;
best_geometry.ts = best_geometry.Rt*best_geometry.t;

best_geometry.b = extract_dimension(0, 'c')/(N+1);
best_geometry.h = best_geometry.b*best_geometry.Rb;
best_geometry.d = best_geometry.Rdh*best_geometry.h;
%}

load('n30_best_geometry')
bg = best_geometry;
figure
scatter(bg.L_distribution(2:end-1), bg.Tr_distribution(1:end-1));
hold on
stairs(bg.L_distribution, [bg.t_distribution(1:end-1), bg.t_distribution(end-1)]);
stairs(bg.L_distribution(1:end), bg.bottom_t);

xlabel('Wing Span (m)')
ylabel('Thickness (mm)')
ylim([0 10]);

legend('Rib Thickness', 'Top Skin Thickness', 'Bottom Skin Thickness')

%scatter(bg.L_distribution(1:end), best_n.t_distribution, 'x');
%legend('Rib Thickness', 'Top Skin Thickness', 'Bottom Skin Thickness')
%}
%{
current_n = T(T.20==i, :);
chosen = T(T.total_weight == min(current_n.total_weight), :)
scatter(chosen.L_distribution(2:end), bg.Tr_distribution)
%}
total_weight = 2700*(bg.total_weight + bg.bottom_weight)/1000

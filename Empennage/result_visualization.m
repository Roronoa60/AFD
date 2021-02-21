close all

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
total_weight = cell2mat(weights(:, 9));
As_bt = cell2mat(weights(:, 10));
te = cell2mat(weights(:,14));
T = table(n, t, Rt, Rdh, Rb, F, L, root_weight, total_weight, As_bt, te);

figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.F, current_n.total_weight);
    hold on
end
legend(string(12:2:34));
xlabel('F');
ylabel('total weight');

%scatter(F, total_weight);

figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.root_weight, current_n.total_weight);
    hold on
end
legend(string(4:2:40));
xlabel('root weight');
ylabel('total weight');

%scatter(root_weight, total_weight);

figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.L, current_n.total_weight);
    hold on
end
legend(string(4:2:40));
xlabel('L');
ylabel('total weight');

figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.n, current_n.total_weight);
    hold on
end

legend(string(4:2:40));
xlabel('n');
ylabel('total weight');

figure
for i=4:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.te, current_n.total_weight);
    hold on
end

legend(string(4:2:40));
xlabel('te');
ylabel('total weight');

figure
for i=[4, 8, 12, 14, 20, 26, 32, 40]
    current_n = T(T.n==i, :);
    scatter(current_n.As_bt, current_n.total_weight);
    hold on
end

legend(string([4, 8, 12, 14, 20, 26, 32, 40]));
xlabel('As_bt');
ylabel('total weight');

sorted_weights = sortrows(T, 'total_weight', 'ascend');

bg = best_geometry;
figure
scatter(bg.L_distribution(2:end), bg.Tr_distribution);
ylabel('Rib Thickness (mm)')
hold on
scatter(bg.L_distribution(1:end), bg.t_distribution);
xlabel('Wing Span (m)')
ylabel('Thickness (mm)')


scatter(bg.L_distribution(1:end), best_n.t_distribution, 'x');

legend('Rib Thickness', 'Top Skin Thickness', 'Bottom Skin Thickness')

total_weight = 2700*(best_n.stringer_weight + bg.total_weight)/1000

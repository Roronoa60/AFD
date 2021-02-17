close all

addpath('Stringer Panel Result Data')
load('varying_L_weights_array.mat')
load('varying_L_best_geometry.mat')
load('bottom_stringer_for_n40')

titles = {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'L', 'root weight', 'total weight'};

weights = weights_array;
n = weights(:, 1);
t = weights(:, 2);
Rt = weights(:, 3);
Rdh = weights(:, 4);
Rb = weights(:, 5);
F = weights(:, 6);
L = weights(:, 7);
root_weight = weights(:, 8);
total_weight = weights(:, 9);
T = table(n, t, Rt, Rdh, Rb, F, L, root_weight, total_weight);

figure
for i=12:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.F, current_n.total_weight);
    hold on
end
legend(string(12:2:34));
xlabel('F');
ylabel('total weight');

%scatter(F, total_weight);

figure
for i=12:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.root_weight, current_n.total_weight);
    hold on
end
legend(string(12:2:40));
xlabel('root weight');
ylabel('total weight');

%scatter(root_weight, total_weight);

figure
for i=12:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.L, current_n.total_weight);
    hold on
end
legend(string(12:2:40));
xlabel('L');
ylabel('total weight');

figure
for i=12:2:40
    current_n = T(T.n==i, :);
    scatter(current_n.n, current_n.total_weight);
    hold on
end

legend(string(12:2:40));
xlabel('n');
ylabel('total weight');
ylim([0 1000])
sorted_weights = sortrows(T, 'total_weight', 'ascend');


bg = best_geometry;
figure
scatter(bg.L_distribution(2:end), bg.Tr_distribution);
ylabel('Rib Thickness (mm)')
hold on
scatter(bg.L_distribution(1:end), bg.t_distribution);
xlabel('Wing Span (m)')
ylabel('Thickness (mm)')
ylim([0 10])

scatter(bg.L_distribution(1:end), best_n.t_distribution, 'x');

legend('Rib Thickness', 'Top Skin Thickness', 'Bottom Skin Thickness')

total_weight = 2700*(best_n.stringer_weight + bg.total_weight)/1000

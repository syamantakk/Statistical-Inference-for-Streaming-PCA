% Data
d_values = [100, 500, 1000, 1500, 2000];

oja_times = [0.0077, 0.0168, 0.0338, 0.0604, 0.0731];
oja_errors = [0.0096, 0.0040, 0.0044, 0.0309, 0.0235];

bootstrap1_times = [0.0134, 0.0283, 0.0484, 0.0503, 0.0966];
bootstrap1_errors = [0.0065, 0.0010, 0.0100, 0.0019, 0.0059];

bootstrap10_times = [0.0448, 0.1724, 0.3194, 0.4650, 0.5888];
bootstrap10_errors = [0.0055, 0.0135, 0.0509, 0.0549, 0.0491];

bootstrap20_times = [0.0911, 0.3201, 0.5883, 0.8517, 1.1071];
bootstrap20_errors = [0.0295, 0.0275, 0.0544, 0.0773, 0.0426];

bootstrap50_times = [0.2027, 0.7858, 1.4295, 2.1626, 2.7354];
bootstrap50_errors = [0.0037, 0.0636, 0.0781, 0.2038, 0.2299];

% Plot with error bars
figure;
hold on;
errorbar(d_values, bootstrap1_times, bootstrap1_errors, '-s', 'LineWidth', 2, 'MarkerSize', 10);
errorbar(d_values, bootstrap10_times, bootstrap10_errors, '-d', 'LineWidth', 2, 'MarkerSize', 10);
errorbar(d_values, bootstrap20_times, bootstrap20_errors, '-^', 'LineWidth', 2, 'MarkerSize', 10);
errorbar(d_values, bootstrap50_times, bootstrap50_errors, '-v', 'LineWidth', 2, 'MarkerSize', 10);
errorbar(d_values, oja_times, oja_errors, '-o', 'LineWidth', 2, 'MarkerSize', 10);

% Labels and Title
xlabel('Dimension d', 'FontSize', 25);
xticks(d_values); % Explicitly set x-axis ticks
ylabel('Time (seconds)', 'FontSize', 25);
title('Computation Time Comparison', 'FontSize', 25);

% Legend
legend({'Bootstrap (b=1)', 'Bootstrap (b=10)', ...
    'Bootstrap (b=20)', 'Bootstrap (b=50)', 'OjaVarEst (Our Algorithm)'}, 'FontSize', 25, 'Location', 'northwest');

% Axis and Grid
set(gca, 'FontSize', 25);
grid on;
hold off;
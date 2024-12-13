% Define method names
methods = {'CMVP', 'GMVP', 'LMVP', 'SC1F', 'SCID', 'MDP', 'Equal Weight'};

% Calculate cumulative returns, risks, mean returns, and Sharpe ratios
cumulative_returns = [cr_x0(end), cr_x1(end), cr_x2(end), cr_x4(end), cr_x5(end), cr_xmdp(end), cr_x6(end)];
risks = [std(out_of_sample_returns0), std(out_of_sample_returns1), std(out_of_sample_returns2), ...
         std(out_of_sample_returns4), std(out_of_sample_returns5), std(out_of_sample_returnsmdp), std(out_of_sample_returns6)];
mean_returns = [mean(out_of_sample_returns0)-1, mean(out_of_sample_returns1)-1, mean(out_of_sample_returns2)-1, ...
                mean(out_of_sample_returns4)-1, mean(out_of_sample_returns5)-1, mean(out_of_sample_returnsmdp)-1, mean(out_of_sample_returns6)-1];
sharpe_ratios = mean_returns ./ risks;

% Create table for cumulative returns, risks, mean returns, and Sharpe ratios
performance_table = table(methods', cumulative_returns', risks', mean_returns', sharpe_ratios', ...
    'VariableNames', {'Method', 'Cumulative_Return', 'Risk', 'Mean_Return', 'Sharpe_Ratio'});


% Short positions (assuming x1asp, x2asp, x4asp, x5asp are defined)
short_positions = [NaN, x1asp, x2asp, x4asp, x5asp, NaN, NaN];
short_positions_table = table(methods', short_positions', 'VariableNames', {'Method', 'Average_Short_Position'});

% Ensure VaR_values and CVaR_values have only 7 elements
VaR_values = VaR_values(1:7);
CVaR_values = CVaR_values(1:7);

% Create table for VaR and CVaR values
VaR_table = table(methods', VaR_values', CVaR_values' - 1, 'VariableNames', {'Method', 'VaR', 'CVaR'});

% Write tables to CSV files
writetable(performance_table, 'Performance_Metrics_w2000.csv');

writetable(short_positions_table, 'Short_Positions_w2000.csv');
writetable(VaR_table, 'VaR_CVaR_Values_w2000.csv');

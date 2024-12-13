clear all;

% Data loading

% csv = csvread('y_wetf.csv',1,1)+1;
% ret=csv(1:469,:);

% csv = csvread('y_W500.csv',1,1)+1;
% ret=csv(1:557,:);

% csv = csvread('y_wcur.csv',1,1)+1;
% ret=csv(1:452,:);

% csv = csvread('FF100.csv',1,1)/100+1;
% ret=csv(140:312,:);

% csv = csvread('ywsp98.csv',1,1)+1;
% ret=csv(81:385,:);

csv = csvread('y_wr2000.csv',1,1)+1;
ret=csv(1:419,:);

cor = corr(ret);
window_size =367;%Rolling window size T
ret0=csv(1:window_size,:);
s0=cov(ret0);
n = size(ret,2);
m= size(ret, 1);
e = ones(n,1);
I = eye(n);
%%
xx0 = [];%save optimal weights for each methods
xx1 = [];
xx4 = [];
xx5 = [];
x22 = [];
xx7=[];
xxmdp = [];
ii=[];

% setting rebalance period(monthly)
rebalance_interval =4; % Rebalancing interval
rebalance_counter = 1; % Initialize counter
out_sample_counter = 1;
% 使用 CVX 优化得到初始权重

opts = optimoptions('quadprog','Display','none');
% % Main loop, rolling window
for t = 1:(m - window_size)  % from t=1
    % Select the return data within the current window
    current_data = ret(1:(t + window_size - 1), :);

    %  Calculate weights at each rebalancing
    if mod(rebalance_counter, rebalance_interval) == 1
        % Compute current S
        S = cov(current_data);
        
        scid = covCor(current_data); %SCID shrinkage
        s_le= cov1Para(current_data);% SC1F shrinkage
       
        % solve GMVP
        cvx_begin quiet
            variable x1(size(ret, 2))
            minimize(0.5 * quad_form(x1, S))
            sum(x1) == 1;
           
        cvx_end
        xx1 = [xx1 x1];

        % solve CMVP
        cvx_begin quiet
            variable x(size(ret, 2))
            minimize(0.5 * quad_form(x, S))
            sum(x) == 1;
            x >= 0;
        cvx_end
        xx0 = [xx0 x];

        % solve SC1F
        cvx_begin quiet;
        variable x4(n)
        minimize(x4' * s_le * x4)
        subject to
            e' * x4 == 1;

        cvx_end
        xx4 = [xx4 x4];

        % solve SCID
        cvx_begin quiet;
        variable x5(n)
        minimize(x5' * scid * x5)
        subject to
            e' * x5 == 1;
            
        cvx_end
        xx5 = [xx5 x5];

        xew = 1/n * ones(n, 1); %Equal-weight

        % solve MDP
        mu = mean(current_data)';
        sigma = std(current_data)';
        MDPprob = optimproblem('ObjectiveSense','minimize');
        nX = size(mu,1);

        y = optimvar('y',nX,1,'LowerBound',0); 
        tau = optimvar('tau',1,1,'LowerBound',0); 

        sigma = std(current_data)'; % 标准差
        MDPprob.Constraints.sigmaSumToOne = sigma'*y == 1;
        MDPprob.Constraints.sumToTau = sum(y) == tau;

        %objective function for MDP
        MDPprob.Objective = y'*S*y;

        % solve
        [wMDP,fMDP] = solve(MDPprob,'options',opts);
        xMDP = wMDP.y/wMDP.tau;

        xxmdp = [xxmdp, xMDP];

    end

    % Calculate out-of-sample returns
    out_of_sample_returns1(out_sample_counter) = ret(t + window_size, :) * x1;
    out_of_sample_returns0(out_sample_counter) = ret(t + window_size, :) * x;
    out_of_sample_returns4(out_sample_counter) = ret(t + window_size, :) * x4;
    out_of_sample_returns5(out_sample_counter) = ret(t + window_size, :) * x5;
    out_of_sample_returns6(out_sample_counter) = ret(t + window_size, :) * xew;
    out_of_sample_returnsmdp(out_sample_counter) = ret(t + window_size, :) * xMDP;
    
    %Increment the rebalancing counter
    rebalance_counter = rebalance_counter + 1;
    
     % Increment the out-of-sample return counter
    out_sample_counter = out_sample_counter + 1;
end
%%
ii = [];
x22 = [];
yy = [];
best_mu_history = [];  % interval for optimal \rho after each iteration
best_cum_return_history = [];  % best_cum_return
out_sample_counter = 1;  
rebalance_counter = 1;  

for t = 1:(m - window_size)  % t
    current_data = ret(t:(t + window_size - 1), :);  % current data
    

    if mod(rebalance_counter, rebalance_interval) == 1
        S = cov(current_data); 
        cvx_begin quiet
            variable x(n)
            dual variables y z;
            minimize(0.5 * quad_form(x, S))
            subject to
                x' * e == 1;  
                y : x >= 0;  % nonnegative constraints
        cvx_end
        
        
        l1 = y;  % save dual variable
        yy = [yy, l1];
        
        %Compute \rho by lemma2 
        r = l1' * (S \ e);  % inverse S
        p = l1' * (S \ l1) * (e' * (S \ e));
        
        % Interval of \rho
        kk2 = (sqrt(p) + r) / (r^2 - p);
        kk3 = (r - sqrt(p)) / (r^2 - p);
       
        % Divide equally
        mu_range = linspace (kk2/1.1, kk3/1.1, 21);
        
        % choosing the optimal \rho
        best_mu = NaN;
        best_cum_return = -Inf;  % 初始化累计收益为负无穷大
        best_x2 = [];

    % Precompute reusable terms outside of the loop
    % Precompute reusable terms outside of the loop
    S_inv = pinv(S);  % Inverse of the covariance matrix S
    S_inv_e = S_inv*e;  % S \ e is constant for all \rho
    e_S_inv_e = e' * S_inv_e;  % (e' * (S \ e)) is constant for all \rho

    % Compute the first constant term: (S_inv * e) / (e' * S_inv * e)
    term_const = S_inv_e / e_S_inv_e;
   

    % Precompute the second term inside the parentheses (without mu):
    term_mu_part = 2 * (S_inv * l1) - (2 * (e' * (S_inv * l1)) / e_S_inv_e) * S_inv_e;
    term_cr1= current_data*term_const;
    term_cr2= current_data*term_mu_part;
    %\boldsymbol{x}_\rho^* = 2\rho \boldsymbol{\Sigma}^{-1} \boldsymbol{l} + \frac{1 - 2 \rho\boldsymbol{1}' \boldsymbol{\Sigma}^{-1}\boldsymbol{l}}{\boldsymbol{1}' \boldsymbol{\Sigma}^{-1} \boldsymbol{1}} \boldsymbol{\Sigma}^{-1} \boldsymbol{1}.
         
    for mu = mu_range
    % Calculate the first part that depends on mu
        x2 = term_const + mu * term_mu_part;
        portfolio_returns = term_cr1 + mu * term_cr2;
        % Calculate portfolio returns
       
        
        
    % Compute cumulative return
     cumulative_return = cumprod(portfolio_returns);
    mean_return= mean(portfolio_returns);
    std_return = std(portfolio_returns);
    sr_return =  mean_return/std_return;
        % Compare cumulative return
        last_day_return = mean(cumulative_return);
      
    
    % Find the mu that gives the maximum cumulative return
         if last_day_return >=best_cum_return
        best_cum_return = last_day_return;
        best_mu = mu;
        best_x2=x2;%adjust_weights(x2);%x2;%x2;%%%
        end
    end
    %      l_c = best_mu*l1;   
    % s_new = S-(e * l_c' + l_c * e');
    % % [v d]=eig(s_new);
    % %     s_new = v'*abs(d)*v;
    %     cvx_begin sdp quiet
    %         variable x2(n)
    %         minimize(0.5 * quad_form(x2, s_new))
    %         subject to
    %             e' * x2 == 1; 
    %             % 权重总和为1
    %             %sum(min(x2, 0)) >= -0.5; % 非负约束
    %     cvx_end
    %     best_x2 =x2;
        x22 = [x22, best_x2];
        ii = [ii, best_mu];  % save optimal \rho

        % save data 
        best_mu_history = [best_mu_history, best_mu];
        best_cum_return_history = [best_cum_return_history, best_cum_return];
    end
    
    
    if (t + window_size) <= m
        out_of_sample_returns2(out_sample_counter) = ret(t + window_size, :) * best_x2;
        out_sample_counter = out_sample_counter + 1;
    end
    
    rebalance_counter = rebalance_counter + 1; 
end


%% out sample return
cr_x1=cumprod(out_of_sample_returns1);
cr_x0=cumprod(out_of_sample_returns0);
cr_x2=cumprod(out_of_sample_returns2);
cr_x4=cumprod(out_of_sample_returns4);
cr_x5=cumprod(out_of_sample_returns5);
cr_x6=cumprod(out_of_sample_returns6);
cr_xmdp = cumprod(out_of_sample_returnsmdp);

%%
%display cumulative return
disp('Cumulative return:');
disp(cr_x0(end));
disp(cr_x1(end));
disp(cr_x2(end));

disp(cr_x4(end));
disp(cr_x5(end));
disp(cr_xmdp(end));
disp(cr_x6(end));


%% Risk
disp('Risk:');
disp(std(out_of_sample_returns0));
disp(std(out_of_sample_returns1));
disp(std(out_of_sample_returns2));
disp(std(out_of_sample_returns4));
disp(std(out_of_sample_returns5));
disp(std(out_of_sample_returnsmdp));
disp(std(out_of_sample_returns6));

%% Expected Return
disp('Expected Return:');
disp(mean(out_of_sample_returns0)-1);
disp(mean(out_of_sample_returns1)-1);
disp(mean(out_of_sample_returns2)-1);
disp(mean(out_of_sample_returns4)-1);
disp(mean(out_of_sample_returns5)-1);
disp(mean(out_of_sample_returnsmdp)-1);
disp(mean(out_of_sample_returns6)-1);

%%
disp(['Sharpe Ratio:']);
disp((mean(out_of_sample_returns0)-1)/std(out_of_sample_returns0));
disp((mean(out_of_sample_returns1)-1)/std(out_of_sample_returns1));
disp((mean(out_of_sample_returns2)-1)/std(out_of_sample_returns2));
disp((mean(out_of_sample_returns4)-1)/std(out_of_sample_returns4));
disp((mean(out_of_sample_returns5)-1)/std(out_of_sample_returns5));
disp((mean(out_of_sample_returnsmdp)-1)/std(out_of_sample_returnsmdp));
disp((mean(out_of_sample_returns6)-1)/std(out_of_sample_returns6));

%%
figure(11)
plot(cr_x0,'g','DisplayName','CMVP');
hold on;
plot(cr_x1,'b','DisplayName','GMVP');
hold on;
plot(cr_x2,'r','DisplayName','LMVP');
hold on;

plot(cr_x4,'y','DisplayName','scvf');
hold on;
plot(cr_x5,'m','DisplayName','scid');
hold on;
plot(cr_xmdp,'r-o','DisplayName','MDP');
hold on;
plot(cr_x6,'b--','DisplayName','ew');


legend({'CMVP','GMVP','LMVP','SC1F','SCID','MDP','EW'},...
    'Location','northwest')

%%
%short selling
x1asp=sum(sum(min(xx1,0)))/size(xx1,2)/n;
x2asp=sum(sum(min(x22,0)))/size(xx1,2)/n;
x4asp=sum(sum(min(xx4,0)))/size(xx1,2)/n;
x5asp=sum(sum(min(xx5,0)))/size(xx1,2)/n;

disp(['Average short position for x1 (GMVP): ', num2str(x1asp)]);
disp(['Average short position for x2 (LMVP): ', num2str(x2asp)]);
disp(['Average short position for x4 (SC1F): ', num2str(x4asp)]);
disp(['Average short position for x5 (SCID): ', num2str(x5asp)]);


%%
% Confidence level for VaR and CVaR
confidence_level = 0.95;  % 95% confidence level


% Define method names
% Preallocate arrays for VaR and CVaR values
VaR_values = zeros(1, 8);
CVaR_values = zeros(1, 8);

% Define method names
methods = {'CMVP', 'GMVP', 'LMVP', 'SC1F', 'SCID', 'MDP', 'Equal Weight','x2_trace'};

% Calculate VaR and CVaR for each method
VaR_values(1) = computeHistoricalVaR(out_of_sample_returns0,0.95,0)-1;
CVaR_values(1) = calculateCVaR(out_of_sample_returns0, confidence_level);

VaR_values(2) = computeHistoricalVaR(out_of_sample_returns1,0.95,0)-1;
CVaR_values(2) = calculateCVaR(out_of_sample_returns1, confidence_level);

VaR_values(3) = computeHistoricalVaR(out_of_sample_returns2,0.95,0)-1;
CVaR_values(3) = calculateCVaR(out_of_sample_returns2, confidence_level);

VaR_values(4) = computeHistoricalVaR(out_of_sample_returns4,0.95,0)-1;
CVaR_values(4) = calculateCVaR(out_of_sample_returns4, confidence_level);

VaR_values(5) = computeHistoricalVaR(out_of_sample_returns5,0.95,0)-1;
CVaR_values(5) = calculateCVaR(out_of_sample_returns5, confidence_level);

VaR_values(6) = computeHistoricalVaR(out_of_sample_returnsmdp,0.95,0)-1;
CVaR_values(6) = calculateCVaR(out_of_sample_returnsmdp, confidence_level);

VaR_values(7) = computeHistoricalVaR(out_of_sample_returns6,0.95,0)-1;
CVaR_values(7) = calculateCVaR(out_of_sample_returns6, confidence_level);


VaR_values(8) = computeHistoricalVaR(out_of_sample_returns6,0.95,0)-1;
CVaR_values(8) = calculateCVaR(out_of_sample_returns6, confidence_level);


% Add a dashed line separator
disp('--------------------------------------');

% Display all VaR values together
disp('VaR values for all methods:');
for i = 1:length(methods)
    disp(['VaR ' methods{i} ': ' num2str(VaR_values(i))]);
end

% Add a dashed line separator
disp('--------------------------------------');

% Display all CVaR values together
disp('CVaR values for all methods:');
for i = 1:length(methods)
    disp(['CVaR ' methods{i} ': ' num2str(CVaR_values(i)-1)]);
end


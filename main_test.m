% ========================================================================
% MATLAB Code for "Modeling and Optimizing Dual-Connectivity Activation in Cellular Networks"
% 
% Copyright (c) 2025 Caglar Tunc, Sabanci University
% All rights reserved.
% 
% This code is provided for research and academic purposes only. 
% Please cite this paper if you use (any part of) this code for your
% research.

run gen_config.m
run define_matrices.m

A1s = cell(1,K);
A2s = cell(1,K);
Ls = cell(1,K);
L0s = cell(1,K);
L1s = cell(1,K);
L2s = cell(1,K);
P = cell(1,K);

for m = 1:K
    Q = Qs{m};
    R = Rs{m};

    P1 = (linsolve([Q ones(S,1)]',[zeros(1,S) 1]'))';
    P(m) = {P1};
    temp = Rs{m};

    offset = 0;

    A = Q/R;
    [Z1,D1] = schur(A);
    cluster = zeros(1,length(D1(1,:)));
    zero = 0; % # of zero eig.values
    neg = 0; % # of negative eig.values
    pos = 0; % # of positive eig.values
    eps = 1e-14;
    for k = 1:size(D1,1)
        if(D1(k,k)<-eps)
            cluster(k) = 2;
            neg = neg+1;
        elseif(D1(k,k)>eps)
            cluster(k) = 1;
            pos = pos+1;
        else %D1(k,k)=0
            cluster(k) = 3;
            zero = zero+1;
        end
    end
    [Z,D] = ordschur(Z1,D1,cluster);

    A_zero = D(1:zero,1:zero);
    A_minus = D(zero+1:zero+neg,zero+1:zero+neg);
    A_plus = D(zero+neg+1:end,zero+neg+1:end);
    A1 = D(1:zero,zero+1:zero+neg);
    A2 = D(1:zero,zero+neg+1:end);
    A3 = D(zero+1:zero+neg,zero+neg+1:end);

    X1 = -D(1:zero,zero+1:end)/D(zero+1:end,zero+1:end);
    X2 = lyap(-A_minus,A_plus,A3);

    M1 = [eye(size(X1,1)) -X1; zeros(size(X1,2),size(X1,1)) eye(size(X1,2))];
    M2 = [eye(S-size(X2,1)-size(X2,2),S-size(X2,1)-size(X2,2)) zeros(S-size(X2,1)-size(X2,2),size(X2,1)) zeros(S-size(X2,1)-size(X2,2),size(X2,2));
        zeros(size(X2,1),S-size(X2,1)-size(X2,2)) eye(size(X2,1)) -X2;
        zeros(size(X2,2),S-size(X2,1)-size(X2,2)) zeros(size(X2,2),size(X2,1)) eye(size(X2,2))];

    Y = Z*M1*M2;
    Y(abs(Y)<1e-14) = 0;
    Y1 = inv(Y);

    YAY = Y1*A*Y;

    A1 = YAY(zero+1:zero+neg, zero+1:zero+neg);
    A2 = YAY(zero+neg+1:end, zero+neg+1:end);

    L0 = Y1(1:zero,:);
    L1 = Y1(zero+1:zero+neg,:);
    L2 = Y1(zero+neg+1:end,:);

    A1s(m) = {A1};
    A2s(m) = {A2};
    Ls(m) = {Y1};
    L0s(m) = {L0};
    L1s(m) = {L1};
    L2s(m) = {L2};

end

% Solving for boundaries

H = [];
v = (2*K+1)*S;

% Conditions 1-4

counter = 0;
cs = ones(K+1,S);
c_count = S*ones(1,K+1);
for k = 1:K+1
    for in = 1:S
        s = S+1-in;
        if(k==1) %c0 conditions (Condition 1)
            if(d(1,s)>0)
                cs(k,s) = 0;
                Q10 = [Q10(1:s-1,:);Q10(s+1:end,:)];
                counter = counter+1;
                c_count(k) = c_count(k) - 1;
            end
        elseif(k==K+1) %cK conditions (Condition 4)
            if(d(K,s)<0)
                cs(k,s) = 0;
                temp = Q1s{k-1};
                temp = [temp(1:s-1,:);temp(s+1:end,:)];
                Q1s(k-1) = {temp};
                counter = counter+1;
                c_count(k) = c_count(k) - 1;
            end
        elseif(k<=K) %ck conditions (Conditions 2&3)
            if(((d(k-1,s)>0 && d(k,s)>0) || (d(k-1,s)<0 && d(k,s)<0)) || ((d(k-1,s)<0 && d(k,s)>0) && (d1(k-1,s)>0 || d1(k-1,s)<0)))
                cs(k,s) = 0;
                temp = Q1s{k-1};
                temp = [temp(1:s-1,:);temp(s+1:end,:)];
                Q1s(k-1) = {temp};
                counter = counter+1;
                c_count(k) = c_count(k) - 1;
            end
        end
    end
end


%-------------------------------------------------------------------------

% Condition 5 - c(0) & a(1)
A1 = A1s{1};
A2 = A2s{1};
L = Ls{1};
L0 = L0s{1};
L1 = L1s{1};
L2 = L2s{1};

count_zeros = S-length(A1)-length(A2);
M1 =  [eye(count_zeros)*L0;expm(-A1*T0)*L1;expm(-A2*T(1))*L2]*Rs{1};
temp = [-Q10;M1;zeros((K-1)*S+sum(c_count(2:end)),S)];
H = [H temp];
%-------------------------------------------------------------------------


for k = 1:K-1
    counter = 0;

    % Condition 6 - a(k) & c(k) & a(k+1)

    A1k = A1s{k};
    A1k1 = A1s{k+1};
    A2k = A2s{k};
    A2k1 = A2s{k+1};

    Lk = Ls{k};
    L0k = L0s{k};
    L1k = L1s{k};
    L2k = L2s{k};
    Lk1 = Ls{k+1};
    L0k1 = L0s{k+1};
    L1k1 = L1s{k+1};
    L2k1 = L2s{k+1};

    Q1 = Q1s{k};
    Rk = Rs{k};
    Rk1 = Rs{k+1};
    pre = counter;
    if(k==1)
        delta1 = T(1)-T0;
    else
        delta1 = T(k)-T(k-1);
    end
    delta2 = T(k+1)-T(k);
    count_zeros1 = S-length(A1k)-length(A2k);
    count_zeros2 = S-length(A1k1)-length(A2k1);

    M1 =  [eye(count_zeros1)*L0k;expm(A1k*delta1)*L1k;expm(0*A2k)*L2k]*Rk;
    M2 =  [eye(count_zeros2)*L0k1;expm(0*A1k1)*L1k1;expm(-A2k1*delta2)*L2k1]*Rk1;

    temp = [zeros((k-1)*S+sum(c_count(1:k)),S); -M1; -Q1; M2; zeros(S*(K-k-1)+sum(c_count(k+2:end)),S)];
    H = [H temp];

end
%-------------------------------------------------------------------------


% Condition 9 - a(K) & c(K)

A1 = A1s{K};
A2 = A2s{K};
L = Ls{K};

count_zeros = S-length(A1)-length(A2);

if(K>1)
    temp = [eye(count_zeros) zeros(count_zeros,length(A1)+length(A2)); zeros(length(A1),count_zeros) expm(A1*(B-T(K-1))) zeros(length(A1),length(A2)); zeros(length(A2),count_zeros+length(A1)) expm(-A2*0)];
else
    temp = [eye(count_zeros) zeros(count_zeros,length(A1)+length(A2)); zeros(length(A1),count_zeros) expm(A1*(B)) zeros(length(A1),length(A2)); zeros(length(A2),count_zeros+length(A1)) expm(-A2*0)];
end

M1 = temp*L*Rs{K};
Q1 = Q1s{K};
temp = [zeros((K-1)*S+sum(c_count(1:end-1)),S);M1;Q1];
H = [H temp];

% Condition 10

temp = zeros(c_count(1),S);
index = find(cs(1,:)==1);
for k = 1:length(index)
    temp(k,index(k)) = 1;
end
for k = 1:K
    A1 = cell2mat(A1s(k));
    A2 = cell2mat(A2s(k));
    L0 = cell2mat(L0s(k));
    L1 = cell2mat(L1s(k));
    L2 = cell2mat(L2s(k));
    count_zeros = S-length(A1)-length(A2);

    if(k==1)
        T1 = T0;
    else
        T1 = T(k-1);
    end
    T2 = T(k);

    r1 = L0*(T2-T1);
    r2 = A1\(expm(A1*(T2-T1))-eye(size(A1)))*L1;
    r3 = expm(-A2*T2)*inv(A2)*(expm(A2*(T2-T1)))*L2;

    M = [r1;r2;r3];

    temp = [temp;M];

    temp1 = zeros(c_count(k+1),S);
    index = find(cs(k+1,:)==1);
    for m = 1:length(index)
        temp1(m,index(m)) = 1;
    end

    temp = [temp;temp1];
end
M1 = temp*ones(S,1);
H = [H M1];
%-------------------------------------------------------------------------

rhs = [zeros(1,size(H,2)-1) 1];
z = rhs/H;

c0 = z(1:c_count(1));
c = cell(1,K);
as = zeros(K,S);
index = c_count(1)+1;
for k = 1:K
    as(k,:) = z(index:index+S-1);
    index = index+S;
    c(k) = {z(index:index+c_count(k+1)-1)};
    index = index+c_count(k+1);
end

c_all0 = zeros(1,S);

index = find(cs(1,:)==1);
temp = c0;
for m = 1:length(index)
    c_all0(index(m)) = temp(m);
end

c_all = zeros(K,S);

for k = 1:K
    index = find(cs(k+1,:)==1);
    temp = c{k};
    for m = 1:length(index)
        c_all(k,index(m)) = temp(m);
    end
end

ignored_boundary_states = [2,3,5,6];
ignored_boundaries = find(T==Td);
normalization_factor = 1;
for bIdx = ignored_boundaries
    for sIdx = ignored_boundary_states
        normalization_factor = normalization_factor - c_all(bIdx,sIdx);
        c_all(bIdx,sIdx) = 0;
    end
end
normalization_factor;

c_all0 = c_all0./normalization_factor;
c_all = c_all./normalization_factor;

% Solve f(x)
f = cell(1,K);
F = cell(1,K);

for k = 1:K
    A1 = A1s{k};
    A2 = A2s{k};
    L = Ls{k};
    count_zeros = S-length(A1)-length(A2);
    if(k==1)
        f(k) = {@(x) as(k,:)*[eye(count_zeros) zeros(count_zeros,length(A1)+length(A2)); zeros(length(A1),count_zeros) expm(A1*x) zeros(length(A1),length(A2)); zeros(length(A2),count_zeros+length(A1)) expm(-A2*(T(k)-x))]*L/normalization_factor};
    else
        f(k) = {@(x) as(k,:)*[eye(count_zeros) zeros(count_zeros,length(A1)+length(A2)); zeros(length(A1),count_zeros) expm(A1*(x-T(k-1))) zeros(length(A1),length(A2)); zeros(length(A2),count_zeros+length(A1)) expm(-A2*(T(k)-x))]*L/normalization_factor};
    end
end

% Calculating PDF and CDF before removing transmission states
delta = config.params.delta;
xs = cell(1,K);
pre = -delta;
for k = 1:K
    xs{k} = (pre+delta):delta:T(k);
    pre = T(k);
end
fs = zeros(S,length(0:delta:B));
Fs = zeros(S,length(0:delta:B));
Fs(:,1) = c_all0';
summ = zeros(1,S);
summ = summ + c_all0;


iter = 1;
Fs(:,iter) = c_all0';
for k = 1:K % going through each regime
    this_pdf = f{k}; % PDF of this regime
    for xIdx = 1:(length(xs{k}))
        x = xs{k}(xIdx);
        pdf_val = max(0,this_pdf(x)); % to avoid very small negative values
        fs(:,iter) = pdf_val;
        if(k==1 && xIdx==1)
            Fs(:,iter) = Fs(:,iter) + fs(:,iter);
        else
            Fs(:,iter) = Fs(:,iter-1) + delta/2*(fs(:,iter)+fs(:,iter-1));
        end
        if(xIdx == length(xs{k}))
            Fs(:,iter) = Fs(:,iter) + c_all(k,:)';
        end
        iter = iter + 1;
    end
end

% Adjusting PDF and CDF after removing transmission states
transmission_states = 7:12;
idle_states = setdiff(1:S,transmission_states);
total_idle_probability = sum(Fs(idle_states,end));
total_tx_probability = sum(Fs(transmission_states,end));

fs_idle = fs(idle_states,:)./total_idle_probability;
Fs_idle = Fs(idle_states,:)./total_idle_probability;

P_m = 1;
P_s = 3;
power_vec = (P_m + P_s)*ones(1,6); % power consumption vector for each state [(1,0), (1,1), (1,2), (2,0), (2,1), (2,2)]
power_vec(1) = P_m;
power_vec(4) = P_m;
avg_power = sum(Fs_idle(:,end)'.*power_vec);
buffer_cdf = sum(Fs_idle,1);
mean_buffer = find(buffer_cdf>=0.5,1)*config.params.delta;

% Calculate SCG activation and deactivation rates
time_act = 0;
time_deact = 0;
if Ta~=0
    rate_act = 0;
    rate_de = 0;

    pi_mat = [K,S]; % steady-state prob. vector for K regimes and S states
    S0 = [1, 4, 7, 10]; % SCG off state indices
    S1 = setdiff(1:S,S0); % SCG on state indices
    r_hat = zeros(1,K);
    r_bar = zeros(1,K);


    time_act = 0;
    time_deact = 0;
    for k = 1:4
        Q = Qs{k};
        R = Rs{k};

        if(T(k)>Td && T(k)<=Ta) % special treatment for regime 2 since the MC is not irreducible
            Q0 = Q(S0,S0);
            P0 = (linsolve([Q0 ones(length(S0),1)]',[zeros(1,length(S0)) 1]'))';
            Q1 = Q(S1,S1);
            P1 = (linsolve([Q1 ones(length(S1),1)]',[zeros(1,length(S1)) 1]'))';
            r_regime = diag(Rs{k});
            r_hat(k) = P0*r_regime(S0);
            r_bar(k) = P1*r_regime(S1);
        else
            P1 = (linsolve([Q ones(S,1)]',[zeros(1,S) 1]'))';
            P1(abs(P1)<1e-14) = 0;
            P1 = P1./sum(P1);
            pi_hat = P1(S0)./sum(P1(S0));
            pi_bar = P1(S1)./sum(P1(S1));
            r_regime = diag(Rs{k});
            r_hat(k) = pi_hat*r_regime(S0);
            r_bar(k) = pi_bar*r_regime(S1);
        end

        this_pdf = f{k};
        l = length(buffer_cdf);
        if T(k) == Td
            buffer_prob = buffer_cdf(round(l/B*Td));
            this_cdf = buffer_cdf(1:round(l/B*Td))./buffer_prob;
            [~,index] = find(this_cdf>=0.5,1);
            this_x = xs{1};
            regime_mean = this_x(index);
            time_act = time_act + (Td-regime_mean)/r_hat(k);
            % elseif k == K
        elseif T(k) > Ta % if above T(a)
            buffer_prob = buffer_cdf(end)-buffer_cdf(round(l/B*Ta));
            this_cdf = (buffer_cdf(round(l/B*Ta):end)-buffer_cdf(round(l/B*Ta)))./buffer_prob;
            [~,index] = find(this_cdf>=0.5,1);
            regime_mean = (B-Ta)/length(this_cdf)*index + Ta;
            time_deact =  time_deact + (regime_mean)/r_bar(k);
        end
    end
    if(Ta ~= Td)
        time_act = (Ta-Td)/r_hat(2) + time_act;
        time_deact = -(Ta-Td)/r_bar(2) - time_deact;
    end
end

% Storing configurations and the output data

if config.save_data
    run gen_save_data.m
    if ~exist(config.output_dir, 'dir')
        % Folder does not exist, so create it
        mkdir(config.output_dir);
        fprintf('Folder "%s" created successfully.\n', config.output_dir);
    end
    save(config.output_dir+config.file_name, 'data');
end

if config.plot_data
    run plot_data.m
end

% Display important statistics

disp(['Avg. activation time: ',num2str(time_act),' ms, avg. deactivation time: ',num2str(time_deact), ' ms'])
disp(['Average delay: ',num2str(mean_buffer/config.params.rate),' ms'])
disp(['Mean power: ',num2str(avg_power)])

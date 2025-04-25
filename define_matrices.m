% Defining Q, and R matrices for K regimes and (K+1) boundaries

K = config.numRegimes;
S = config.numStates;

am1 = config.params.am1;
rm1 = config.params.rm1;
as1 = config.params.as1;
rs1 = config.params.rs1;
ps1 = rs1/(rs1+as1); % probability of state-1 for SCG
ps2 = 1-ps1; % probability of state-2 for SCG

mu_m1 = config.params.mu_m1;
mu_m2 = config.params.mu_m2;
mu_s1 = config.params.mu_s1;
mu_s2 = config.params.mu_s2;

packet_drift = config.params.packet_drift;
ba = config.params.ba;
ta = config.params.ta;

T0 = config.params.T0;
B = config.params.B;
Td = config.params.Td;
Ta = config.params.Ta;
T = config.params.T(2:end);

Qs = cell(1,K);
Q1s = cell(1,K);

        for k = 1:K % Generating Q matrices
            if(T(k)<=Td) % We are in an MCG-only regime
%        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Qs(k) = {[  0       0       0      am1      0        0      ta      0       0       0       0       0; % (1,0,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,1,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,2,0)
           rm1      0       0       0       0        0       0      0       0      ta       0       0; % (2,0,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,1,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,2,0)
           ba       0       0       0       0        0       0      0       0       0       0       0; % (1,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0      ba       0        0       0      0       0       0       0       0; % (2,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0]};%(2,2,1)
            elseif(T(k)<=Ta) % We are in an MCG-only OR MCG+SCG regime
%        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Qs(k) = {[  0       0       0      am1      0        0      ta      0       0       0       0       0; % (1,0,0)
            0       0      as1      0      am1       0       0     ta       0       0       0       0; % (1,1,0)
            0      rs1      0       0       0       am1      0      0      ta       0       0       0; % (1,2,0)
            rm1     0       0       0       0        0       0      0       0      ta       0       0; % (2,0,0)
            0      rm1      0       0       0       as1      0      0       0       0      ta       0; % (2,1,0)
            0       0      rm1      0      rs1       0       0      0       0       0       0      ta; % (2,2,0)
            ba      0       0       0       0        0       0      0       0       0       0       0; % (1,0,1)
            0      ba       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0      ba       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0      ba       0        0       0      0       0       0       0       0; % (2,0,1)
            0       0       0       0      ba        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0       ba       0      0       0       0       0       0]}; % (2,2,1)
            else % We are in an MCG+SCG regime
%        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Qs(k) = {[  0       0       0       0       0        0       0      0       0       0       0       0; % (1,0,0)
            0       0      as1      0      am1       0       0     ta       0       0       0       0; % (1,1,0)
            0      rs1      0       0       0       am1      0      0      ta       0       0       0; % (1,2,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,0,0)
            0      rm1      0       0       0       as1      0      0       0       0      ta       0; % (2,1,0)
            0       0      rm1      0      rs1       0       0      0       0       0       0      ta; % (2,2,0)
            0    ps1*ba   ps2*ba    0       0        0       0      0       0       0       0       0; % (1,0,1)
            0      ba       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0      ba       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0       0    ps1*ba   ps2*ba     0      0       0       0       0       0; % (2,0,1)
            0       0       0       0      ba        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0       ba       0      0       0       0       0       0]};%(2,2,1)
            end
        end

        for k = 1:K % Generating Q1 matrices
            if(T(k)==Td) % We are at the SCG deactivation boundary
%        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Q1s(k) = {[ 0       0       0      am1      0        0      ta      0       0       0       0       0; % (1,0,0)
           10       0       0       0       0        0       0      0       0       0       0       0; % (1,1,0)
           10       0       0       0       0        0       0      0       0       0       0       0; % (1,2,0)
            rm1     0       0       0       0        0       0      0       0      ta       0       0; % (2,0,0)
            0       0       0      10       0        0       0      0       0       0       0       0; % (2,1,0)
            0       0       0      10       0        0       0      0       0       0       0       0; % (2,2,0)
            ba      0       0       0       0        0       0      0       0       0       0       0; % (1,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0      ba       0        0       0      0       0       0       0       0; % (2,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0]};%(2,2,1)
            elseif(T(k)==Ta) % We are at the SCG activation boundary
%        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Q1s(k) = {[ 0       0       0      am1      0        0      ta      0       0       0       0       0; % (1,0,0)
            0       0      as1      0      am1       0       0     ta       0       0       0       0; % (1,1,0)
            0      rs1      0       0       0       am1      0      0      ta       0       0       0; % (1,2,0)
            rm1     0       0       0       0        0       0      0       0      ta       0       0; % (2,0,0)
            0      rm1      0       0       0       as1      0      0       0       0      ta       0; % (2,1,0)
            0       0      rm1      0      rs1       0       0      0       0       0       0      ta; % (2,2,0)
            ba      0       0       0       0        0       0      0       0       0       0       0; % (1,0,1)
            0      ba       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0      ba       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0      ba       0        0       0      0       0       0       0       0; % (2,0,1)
            0       0       0       0      ba        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0       ba       0      0       0       0       0       0]}; % (2,2,1)
            elseif(k==K) % We are at B
%        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Q1s(k) = {[ 0       0       0       0       0        0       0      0       0       0       0       0; % (1,0,0)
            0       0      as1      0      am1       0       0     ta       0       0       0       0; % (1,1,0)
            0      rs1      0       0       0       am1      0      0      ta       0       0       0; % (1,2,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,0,0)
            0      rm1      0       0       0       as1      0      0       0       0      ta       0; % (2,1,0)
            0       0      rm1      0      rs1       0       0      0       0       0       0      ta; % (2,2,0)
            0    ps1*ba   ps2*ba    0       0        0       0      0       0       0       0       0; % (1,0,1)
            0      ba       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0      ba       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0       0    ps1*ba   ps2*ba     0      0       0       0       0       0; % (2,0,1)
            0       0       0       0      ba        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0       ba       0      0       0       0       0       0]};%(2,2,1)
            else % We are at a regular boundary, set the generator matrix to that of the regime
                Q1s(k) = Qs(k);
            end
        end

% Generating Q10, regime 1 lower boundary matrix
if(Td == 0) % we are at the SCG deactivation boundary
    %        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Q10 =     [ 0       0       0      am1      0        0      ta      0       0       0       0       0; % (1,0,0)
           10       0       0       0       0        0       0      0       0       0       0       0; % (1,1,0)
           10       0       0       0       0        0       0      0       0       0       0       0; % (1,2,0)
            rm1     0       0       0       0        0       0      0       0      ta       0       0; % (2,0,0)
            0       0       0      10       0        0       0      0       0       0       0       0; % (2,1,0)
            0       0       0      10       0        0       0      0       0       0       0       0; % (2,2,0)
            ba      0       0       0       0        0       0      0       0       0       0       0; % (1,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0      ba       0        0       0      0       0       0       0       0; % (2,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0];
else
        %        (1,0,0) (1,1,0) (1,2,0) (2,0,0) (2,1,0) (2,2,0) (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
Q10 =    [  0       0       0      am1      0        0      ta      0       0       0       0       0; % (1,0,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,1,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,2,0)
           rm1      0       0       0       0        0       0      0       0      ta       0       0; % (2,0,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,1,0)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,2,0)
           ba       0       0       0       0        0       0      0       0       0       0       0; % (1,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (1,2,1)
            0       0       0      ba       0        0       0      0       0       0       0       0; % (2,0,1)
            0       0       0       0       0        0       0      0       0       0       0       0; % (2,1,1)
            0       0       0       0       0        0       0      0       0       0       0       0]; % (2,2,1)
end
        d = zeros(K,S);
        for k = 1:K % Generating the drift vectors at each space
            if(T(k)<=Td) % We are in an MCG-only regime
                     %    (1,0,0)      (1,1,0)       (1,2,0)      (2,0,0)     (2,1,0)        (2,2,0)   (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
                d(k,:) = [(-mu_m1)        1              1       (-mu_m2)        1              1         packet_drift*ones(1,6)];
            else
                d(k,:) = [(-mu_m1) (-mu_m1-mu_s1) (-mu_m1-mu_s2) (-mu_m2) (-mu_m2-mu_s1) (-mu_m2-mu_s2)   packet_drift*ones(1,6)];
            end
        end

        d1 = zeros(K,S);
        for k = 1:K % Generating the drift vectors at each space
            if(T(k)==Td) % We are at the SCG deactivation boundary
                     %    (1,0,0)      (1,1,0)       (1,2,0)      (2,0,0)     (2,1,0)        (2,2,0)   (1,0,1) (1,1,1) (1,2,1) (2,0,1) (2,1,1) (2,2,1)
                d1(k,:) = [(-mu_m1)        0             0        (-mu_m2)        0               0        packet_drift*ones(1,6)];
            elseif(T(k)==Ta) % We are at the SCG activation boundary
                d1(k,:) = [(-mu_m1) (-mu_m1-mu_s1) (-mu_m1-mu_s2) (-mu_m2) (-mu_m2-mu_s1) (-mu_m2-mu_s2)   packet_drift*ones(1,6)];
            elseif(k==K) % We are at B
                d1(k,:) = [0     (-mu_m1-mu_s1) (-mu_m1-mu_s2)     0    (-mu_m2-mu_s1) (-mu_m2-mu_s2)   0 packet_drift packet_drift 0 packet_drift packet_drift];
            else % We are at a regular boundary, set the generator matrix to that of the regime
                d1(k,:) = d(k,:);
            end
        end

Q10 = Q10 - eye(size(Q10)).*diag(sum(Q10,2));
for k = 1:K
    Qs(k) = {Qs{k} - eye(size(Qs{k})).*diag(sum(Qs{k},2))};
    Q1s(k) = {Q1s{k} - eye(size(Q1s{k})).*diag(sum(Q1s{k},2))};
end
Rs = cell(1,K);
R1s = cell(1,K);
for k = 1:K
    Rs(k) = {diag(d(k,:))};
    R1s(k) = {diag(d1(k,:))};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulated annealing for graph isomorphism %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%%%% Load dataset
load 'GI_dataset.mat'
load 'GI_dataset2.mat'

%% Number of trials
run_c = 10;

%% Problem size
N = 20;

%%%%% Annealing parameters
% T_ini: Initial temperature, T_end: Final temperature, Mcycle: Total cycle number
% tau: Number of cycles before temperature reduction, Cr: Temperature decrease rate
% pn: Node inversion probability
T_ini = 1000;   
T_end = 0.1; 
Mcycle = 40000;
tau = 1;

%% TTS parameters
Ps = 0.99;

%% Data acquisition
aT = zeros(1, run_c);           % Simulation time
correct = 0;                    % Number of correct answers
TTS = 0;  

%%%G1E set
if(N == 5) G1E = G5;
elseif(N == 10) G1E = G10;
elseif(N == 15) G1E = G15;
elseif(N == 20) G1E = G20;
elseif(N == 25) G1E = G25;
elseif(N == 30) G1E = G30;
elseif(N == 35) G1E = G35;
elseif(N == 40) G1E = G40;
elseif(N == 45) G1E = G45;
elseif(N == 50) G1E = G50;
elseif(N == 55) G1E = G55;
elseif(N == 60) G1E = G60;
elseif(N == 65) G1E = G65;
elseif(N == 70) G1E = G70;
elseif(N == 75) G1E = G75;
elseif(N == 80) G1E = G80;
elseif(N == 85) G1E = G85;
elseif(N == 90) G1E = G90;
elseif(N == 95) G1E = G95;
elseif(N == 100) G1E = G100;
end

% G2 which is identical to G1
G2E = G1E;

%%%%% Generating and showing the two graphs
G1 = graph(G1E~=0);
G2 = graph(G2E~=0);

%%%%%% Generation of a Hamiltonian for G1 and G2 in QUBO format
H1 = zeros(N^2,N^2);

    %%%%% Setting penalty parameters
    C1 = 1;
    C2 = 1;

    %%%%%% Adding C1 penalty to the Hamiltonian
    for u=1:N
        for i=1:N
            for v=u:N
                for j=i:N
                    if i==j && u==v
                        H1((u-1)*N+i,(v-1)*N+j) = -C1; 
                    elseif i==j || u==v
                        H1((u-1)*N+i,(v-1)*N+j) = H1((u-1)*N+i,(v-1)*N+j) + C1;
                    end
                end
            end
        end
    end

    %%%%% Adding C2 penalty to the hamiltonian
    for u=1:N
        for v=u:N
            if u~=v
                for i=1:N
                    for j=1:N
                        if i~=j
                            if G2E(u,v) ~= G1E(i,j)
                                H1((u-1)*N+i,(v-1)*N+j) = H1((u-1)*N+i,(v-1)*N+j) + C2;
                            end
                        end
                    end
                end
            end
        end
    end

len = length(H1);
H = H1;

%%%%%% Convertion from QUBO format to ising model for Hamiltonian
J = H/4;

for i=1:len
    for j=1:len
        if i>j
            J(i,j) = J(j,i);
        elseif i==j
            J(i,j) = 0;
        end
    end
end

for i=1:len
     h(i,1) = H(i,i)/2 + sum(J(i,:));
end

h = -h;
J = -J;

%%%%% Find minimum energy
for i = 1:N
    for j = 1:N
        if i == j
            m((i-1)*N+j,1) = 1;
        else
            m((i-1)*N+j,1) = -1;
        end
    end
end
Jm_temp_min = J*m;  
hm_temp_min = transpose(h)*m;
true_min_energy = -sum(Jm_temp_min.*m)/2 - hm_temp_min;

fprintf('SA：N = %d, ',N);
fprintf('EC = %d\n',Mcycle);
for run = 1:run_c

    Cr = power(T_end/T_ini, 1/Mcycle); %%一般的に Cr = 0.8 ~ 0.99
    %%pn = 0.1;

    %%%%% Initializing variables
    % T：temp
    T = zeros(1,Mcycle+1);
    energy = zeros(1,Mcycle+1);
    mi = 2*randi([0,1],len,1)-1;
    mo = zeros(len,Mcycle+1);
    mo(:,1) = mi;

    %%%%% Energy calculation
    Jm_temp = J*mi;  
    hm_temp = transpose(h)*mi;
    energy(1) = -sum(Jm_temp.*mi)/2 - hm_temp;

    %%%%% Start time measurement
    tStart = tic;

    %%%%% Annealing
    for t = 2:Mcycle+1
        
        %%%%% Change temperature
        if t == 1
            T(t) = T_ini;
        elseif rem(t,tau) == 0
            if T_end < T(t-1)
                T(t) = T(t-1) * Cr;
            else
                T(t) = T_ini;
                mi = 2*randi([0,1],len,1)-1;
            end
        else
            T(t) = T(t-1);
        end

        %%%%% Flip a randomly selected node
        select = randi([1 len]);
        mo(:,t) = mi;
        mo(select,t) = -mo(select,t);

        %%%%% Energy calculation
        Jm_temp = J*mo(:,t);  
        hm_temp = transpose(h)*mo(:,t);
        energy(t) = -sum(Jm_temp.*mo(:,t))/2 - hm_temp;

        deltaE = energy(t) - energy(t-1);
        P = exp(-deltaE / T(t));

        if P > rand
            mi(select) = -mi(select);
        else
            energy(t) = energy(t-1);
            mo(select,t) = mi(select);
        end   
    end
    
    %%%%% End time measurement
    aT(1,run) = toc(tStart);

    if min(energy) == true_min_energy
        correct = correct + 1;
        %fprintf('try%d　〇\n',run);
        fprintf('try%d　〇   %.5f[sec]\n',run,aT(1,run));
    else
        %fprintf('try%d　×\n',run);
        fprintf('try%d　×   %.5f[sec]\n',run,aT(1,run));
    end
    
end

%% End of run_c for loop
time = sum(aT)/run_c;  %% Average simulation time per run
acc = correct/run_c;   %% Accuracy rate
if acc ~= 0
    TTS = time * (log(1-Ps) / log(1 - acc));
end
fprintf('Result for  N = %d: Accuracy = %.2f％\n',N,acc*100);
fprintf('Mean time = %.5f[sec], TTS = %.5f[sec]\n',time,TTS);

 


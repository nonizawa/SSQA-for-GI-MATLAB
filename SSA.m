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

%% Annealing parameters
tau = 10;
I0_ini = 1; %%for stochastic
I0_max = 16;
beta = 0.5;
nrnd = 1;
Mcycle = 40000;

%% Data acquisition
aT = zeros(1, run_c);           % Simulation time
correct = 0;                    % Number of correct answers
TTS = 0;                        % Time to solution

%% TTS parameters
Ps = 0.99;

%% Minimum energy
min_energy = 0;

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

fprintf('SSA：N = %d, ',N);
fprintf('EC = %d\n',Mcycle);
for run = 1:run_c
%%%%% Setting parameters for annealing
original = 0;   % 0:Stochasitc, 1:Original(p-bit)
vector = 0;     % Option for original = 0 | 0: serial, 1: vector)

Prnd = 1;

NI = 1 + double(int64(log10(I0_ini/I0_max)/log10(beta)));

%%%%%% Initializing variables
I0 = zeros(1,Mcycle);
energy = zeros(1,Mcycle);

mi = 2*randi([0,1],len,1)-1;
mo = zeros(len,Mcycle+1);
I = zeros(len,Mcycle+1);

vones = ones(len,1);

%%%%% Generation of random signals

rand_store = nrnd*(2*(randi([0,1],len,Mcycle))-1);
rand_store = rand_store.*(rand(len,Mcycle) < Prnd);

%%%%% Compute minimum energy
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
    
%%%%% Start time measurement
tStart = tic;

%%%%%%% Annealing process
for t=1:Mcycle
    
    %%%%% Control of psuedo inverse temperature
    if t==1
        I0(1,t) = I0_ini;
    elseif rem(t,NI*tau)==0  %%I0_iniからI0_maxまでいった時のリセット
        I0(1,t) = I0_ini;
        %%mi = 2*randi([0,1],len,1)-1;
        %%I(:,t) = zeros(len,1);
    elseif rem(t,tau)==0 
        I0(1,t) = I0(1,t-1) / beta;
    else
        I0(1,t) = I0(1,t-1);
    end
    
    %%%%% Energy calculation
    Jm_temp = J*mi;  
    hm_temp = transpose(h)*mi;
    energy(1,t) = -sum(Jm_temp.*mi)/2 - hm_temp;
    
   %%%%%%%%%%%%%%%% Node computation %%%%%%%%%%%%%%%%%%%%%%
   if original == 1
            I(:,t+1) = tanh(I0(1,t) * (h + Jm_temp)) + rand_store(:,t);  %%I0が大きくなることで乱数randに左右されなくなっていく
            mo(:,t+1) = 2*( I(:,t+1)>=0 )-1;
   else
     if vector == 1
        I(:,t+1) =  I(:,t) + h + Jm_temp + rand_store(:,t);
        vI0 = I0(1,t)*vones;
        I(:,t+1) = (I(:,t+1) >= vI0).*(vI0-vones) ...
            + (I(:,t+1) < -vI0).*(-vI0) ...
            + (I(:,t+1) < vI0).*(I(:,t+1) >= -vI0).*I(:,t+1);
        mo(:,t+1) = 2*(I(:,t+1)>=0)-1;
     else
        for a=1:len
        I(a,t+1) =  I(a,t) + h(a,1) + Jm_temp(a,1) + rand_store(a,t);
        if I(a,t+1) >= I0(1,t)
            I(a,t+1) = I0(1,t);
        elseif I(a,t+1) < -I0(1,t)
            I(a,t+1) = -I0(1,t);
        end
        mo(a,t+1) = 2*(I(a,t+1)>=0)-1;
        end
     end
   end
   
   mi = mo(:,t+1); % Transfering node outputs to node inputs in the next cycle
    
end

%%%%% End time measurement
aT(1,run) = toc(tStart);


min_energy = min(min(energy));
if true_min_energy == min_energy
    fprintf('try%d　〇   %.5f[sec]\n',run,aT(1,run));
    correct = correct + 1;
else
    fprintf('try%d　×   %.5f[sec]\n',run,aT(1,run));
end

end %% End of run_c for loop
time = sum(aT)/run_c;  %% Average simulation time per run
acc = correct/run_c;   %% Accuracy rate
TTS = time * (log(1-Ps) ./ log(1 - acc));
fprintf('Result for N = %d: Accuracy = %.2f％\n',N,acc*100);
fprintf('Mean time = %.5f[sec], TTS = %.5f[sec]\n',time,TTS);

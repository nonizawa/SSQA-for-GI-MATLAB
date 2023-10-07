%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Stochastic simulated Quantum annealing for graph isomorphism %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%%%% Load dataset
load 'GI_dataset.mat'
load 'GI_dataset2.mat'

%% Number of trials
run_c = 100;

%% Problem size
N = 20;

%% Annealing parameters
I0 = 2;                         % Inverse temperature 
M = 25;                         % Trotter slices                     
Q_max = 0.5;                    % Maximum inter-layer interaction (minimum: Q_min = 0)
tau = 100;                      % Interval to increase Q (number of cycles)
T = 1;                          % Interference timing
nrnd = 1;                       % Noise amplitude
Prnd = 1;                       % Probability of adding noise
iteration = 4 * tau;            % Number of cycles per iteration
Ni = 4;                        % Number of iterations
Mcycle = M * iteration * Ni;    % Total number of simulation cycles

%% Data acquisition
aT = zeros(1, run_c);           % Simulation time
correct = 0;                    % Number of correct answers
TTS = 0;                        % Time to solution

%% TTS parameters
Ps = 0.99;

%% Minimum energy
min_energy = 0;

%%% G1E set

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


% Generate Graph2, which is identical to Graph1
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

h = zeros(len,1);
for i=1:len
     h(i,1) = H(i,i)/2 + sum(J(i,:));
end

h = -h;
J = -J;

fprintf('SSQA：R = %d, N = %d, ',M,N);
fprintf('EC = %d, SC = %d\n',Mcycle,Mcycle/M);
for run = 1:run_c
%%%%% Setting parameters for annealing
original = 0;   % 0:Stochasitc, 1:Original(p-bit)
vector = 0;     % Option for original = 0 | 0: serial, 1: vector)

%%%%%% Initializing variables
Q = zeros(1,floor(Mcycle/M)+1+1);
energy = zeros(1,floor(Mcycle/M)+1+1,M);
mi = 2*randi([0,1],len,1,M,floor(Mcycle/M)+1+1)-1;
mo = zeros(len,floor(Mcycle/M)+1+1,M);
I = zeros(len,floor(Mcycle/M)+1+1,M);

vones = ones(len,1);

%%%%% Generation of random signals
rand_store = nrnd*(2*(randi([0,1],len,floor(Mcycle/M)+1+1,M))-1);
rand_store = rand_store.*(rand(len,floor(Mcycle/M)+1+1,M) < Prnd);

%%%%% Compute minimum energy
m = zeros(len,1);
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
t = 1;
k = 1;
while t * M <= Mcycle
    for k = 1:M
        if t==1
            Q(1,t) = 0;
        elseif rem(t,tau)==0
            if Q(1,t-1) >= Q_max - 0.0001
                Q(1,t) = 0; %初期化(次のイタレーションへ) 
                if rem(t,2*iteration) == 0 %2Iteration毎に初期化
                    mi(:,1,k,t) = 2*randi([0,1],len,1)-1;
                    I(:,t,k) = zeros(len,1);
                end
            else
            Q(1,t) = Q(1,t-1) + Q_max/((iteration/tau)-1); 
            end
        else
            Q(1,t) = Q(1,t-1);  %%そのまま
        end

    %%%%% Energy calculation
    Jm_temp = J*mi(:,:,k,t);  
    hm_temp = transpose(h)*mi(:,:,k,t);
        if t-T <= 1
            Qm_temp = zeros(len,1);
        elseif k == 1 
            Qm_temp = Q(1,t) * mi(:,:,M,t-T);
        else
            Qm_temp = Q(1,t) * mi(:,:,k-1,t-T);
        end
    energy(1,t,k) = -sum(Jm_temp.*mi(:,:,k,t))/2 - hm_temp;
    
   %%%%%%%%%%%%%%%% Node computation %%%%%%%%%%%%%%%%%%%%%%
   if original == 1
            I(:,t+1,k) = tanh(I0 * (h + Jm_temp + Qm_temp)) + rand_store(:,t,k); %%tanh関数は-1～1  randも-1～１
            mo(:,t+1,k) = 2*( I(:,t+1,k)>=0 )-1;
   else
     if vector == 1
        I(:,t+1,k) =  I(:,t,k) + h + Jm_temp + Qm_temp + rand_store(:,t,k);
        vI0 = I0*vones;
        I(:,t+1,k) = (I(:,t+1,k) >= vI0).*(vI0-vones) ...
            + (I(:,t+1,k) < -vI0).*(-vI0) ...
            + (I(:,t+1,k) < vI0).*(I(:,t+1,k) >= -vI0).*I(:,t+1,k);
        mo(:,t+1,k) = 2*(I(:,t+1,k)>=0)-1;
     else
        for a=1:len
        I(a,t+1,k) =  I(a,t,k) + h(a,1) + Jm_temp(a,1) + Qm_temp(a,1) + rand_store(a,t,k);
        if I(a,t+1,k) >= I0
            I(a,t+1,k) = I0-1;
        elseif I(a,t+1,k) < -I0
            I(a,t+1,k) = -I0;
        end
        mo(a,t+1,k) = 2*(I(a,t+1,k)>=0)-1;
        end
     end
   end
   
   for i = 1:len
       mi(i,1,k,t+1) = mo(i,t+1,k); % Transfering node outputs to node inputs in the next cycle
   end
   
    end
    t = t + 1;
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
fprintf('Result for M = %d, N = %d: Accuracy = %.2f％\n',M,N,acc*100);
fprintf('Mean time = %.5f[sec], TTS = %.5f[sec]\n',time,TTS);

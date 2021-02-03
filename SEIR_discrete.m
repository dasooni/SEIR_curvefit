clear;

%data containing new daily cases in ostergotland, norrkoping
load('seok.mat'); %new daily cases, sept-okt period
load('okno.mat'); %new daily cases, oct-nov period

sn = [seok; okno]; %combined data
avg_sn = movmean(sn,[3 3]); %moving mean for smoothing curve , calculated for 3 days before and 3 days ahead

%initial guess for infection intensity
beta = 0.54;

%start and end point to fit shorter curves, by reducing the period e.g. a
%week
start = 74;
est = 88;

%the objective is to minimize the error between actual data and model
%output
obj = @(b) norm(seir_model2(b)'-avg_sn(start:est));
%solution provided by fminsearch
sol = fminsearch(obj,0.3);
%the difference between model output and actual data
err = obj(sol);
%model output when curve-fit produces a result
res = seir_model2(sol);

%draw model output
plot(res);
hold on
%draw actual data
plot(avg_sn(start:est));


%discrete SEIR-model function
function NC = seir_model2(beta)
    gamma = 1/5; 
    N = 4.62*10^5;
    I0 = 95.5714; 
    E0 = 520;
    T = 15;
    sigma = 1/(5.1);

    S = zeros(1,T);
    E = zeros(1,T);
    I = zeros(1,T);
    R = zeros(1,T);
    NC = zeros(1,T);
    I(1) = I0;
    E(1) = E0;
    S(1) = N-I0-E0;
    NC(1) = I0;
    for t = 1:T-1
    
        dS = (-beta*I(t)*S(t))/N; %suspectible
        dE = (beta*I(t)*S(t))/N - sigma*E(t);
        %dI = (beta*I(t)*S(t)/N - gamma*I(t)); %infected SIR
        dI = sigma*E(t) - (gamma)*I(t); %infected 
        dR = (gamma*I(t));  %removed
    
        %update for each time step
        S(t+1) = S(t) + dS;
        E(t+1) = E(t) + dE;
        I(t+1) = I(t) + dI;
        R(t+1) = R(t) + dR;
        NC(t+1) = dI + dR;
    end
    end

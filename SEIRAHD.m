clear all
close all
clc
% Parameters
%load FB2404
%A=ones(2404)-eye(2404);
alpha = 1/5; %Given
beta = 0.3; %Taken from SEIRDLockdown_NewVersion.m
gamma = 1/10; %(alpha/0.66) = (gamma/0.33)
delta = 1/10; %Given
theta = 1/7; %Given
h = 1/90; %(delta/0.9) = (h/0.1)
omega = 1/63; %(theta/0.9) = (omega/0.1)
Seeds = 100;

%Initialise variables
Asiz = size(A,1);

E = zeros(Asiz,1); %assigns initial seeds
for seed = 1:Seeds
    E(ceil(Asiz*rand),1) = 1;
end

S = ones(Asiz,1) - E; %assigns susceptible population
[I,As,R,H,D] = deal(zeros(Asiz,1)); %all other populations are 0 (no vaccines)
[NewS,NewE,NewI,NewA,NewR,NewH,NewD] = deal(zeros(Asiz,1));

%Mask wearers
Masks = zeros(Asiz,1);
for mask=1:floor(Asiz*0.5) %mask wearing percentage
    Masks(ceil(Asiz*rand)) = 1; %assigns random people to be mask wearers according to percentage defined in for loop
end

%How a Lockdown graph would look like
ALockdown = A;
for i = 1:Asiz 
    for j = i+1:Asiz
        if (rand<0.5) ==1
            ALockdown(i,j) = 0;
            ALockdown(j,i) = 0;
        end
    end
end

AHalfLockdown = A;
for i = 1:Asiz 
    for j = i+1:Asiz
        if (rand<0.25) == 1
            ALockdown(i,j) = 0;
            ALockdown(j,i) = 0;
        end
    end
end

t = 1;
ACurrent = A;
[mask,lockdown] = deal(zeros(200,1));
%Iterations of infection
while (sum(E) + sum(As) + sum(I) > 0) 
    %Population headcounts
    t;
    SumS(t) = sum(S);
    SumE(t) = sum(E);
    SumA(t) = sum(As);
    SumI(t) = sum(I);
    SumH(t) = sum(H);
    if SumH(t) < 0
        SumH(t) = 0;
    end
    SumR(t) = sum(R);
    SumD(t) = sum(D);

    %Lockdown Strategy
    if sum(lockdown) > 60
        ACurrent = AHalfLockdown;
        lockdown(t) = 0.5;
    elseif SumI(t)/Asiz > 0.02 %Upper threshold for lockdown
        ACurrent = ALockdown;
        lockdown(t) = 1;
    elseif (SumI(t)/Asiz) < 0.0005 %Lower threshold for lockdown
        ACurrent = A;
        lockdown(t) = 0;
    else
        lockdown(t) = lockdown(t-1);
    end
    
    
    %Mask Strategy
    if (SumI(t)/Asiz) > 0.001
        mask(t) = 1;
    elseif (SumI(t)/Asiz) < 0.0001
        mask(t) = 0;
    else
        mask(t) = mask(t-1);
    end

    %Mask wearing neighbors
    INeighbors = ACurrent*(or(I,As)); %vector of infected neighbors
    INeighbors_Mask = ACurrent*(or(I,As)).*Masks;
    INeighbors_NMask = INeighbors - INeighbors_Mask;
    
    if mask(t) ==1
        NewE=rand(Asiz,1)<1-(1-(1-0.3*Masks)*beta).^INeighbors_NMask.*(1-(0.05-0.035*Masks)*beta).^INeighbors_Mask; %infected given that at you/neighbor wore mask
    else
        NewE = rand(Asiz,1) < 1 - (1-beta).^(INeighbors);   %infected given that no one wore mask
    end
    %S to E
    NewE = and(NewE,boolean(S)); %people that become E given that they were in S

    %E to I or A
    ERandom = rand(Asiz,1);
    NewI = and((ERandom<alpha), boolean(E));
    NewA = and((ERandom>1-gamma), boolean(E));

    %I to H
    IRandom = rand(Asiz,1);
    NewH = and((IRandom<0.1), boolean(I));
    
    %H to D
    HRandom = rand(Asiz,1);
    anotherrand = rand(Asiz, 1);
    age85 = and ((HRandom < 0.02), boolean(H));
    age85 = and ((anotherrand < 0.3), age85);
    age75 = and ((HRandom > 0.02), (HRandom < 0.06));
    age75 = and (age75, boolean(H));
    age75 = and ((anotherrand < 0.2), age75);
    age65 = and ((HRandom > 0.06), (HRandom < 0.16));
    age65 = and (age65, boolean(H));
    age65 = and ((anotherrand < 0.1), age65);
    ageother = and ((HRandom > 0.16), boolean(H));
    ageother = and ((anotherrand < 0.01), ageother);
    
    NewD = or (age85, age75);
    NewD = or (NewD, age65);
    NewD = or (NewD, ageother);


    %I or H or A to R
    %R had to differently programmed b/c it has multiple in edges
    ItoR = and((IRandom>1-delta), boolean(I));
    HtoR = and(and((HRandom<theta), boolean(H)), not(NewD));
    AtoR = and((rand(Asiz,1)<delta), boolean(As));
    NewR = or(or(ItoR, HtoR), AtoR);

    %Update indicators
    S = S - NewE;
    E = E + NewE - (NewI + NewA);
    I = I + NewI - ItoR - NewH;
    As = As + NewA - AtoR;
    H = H + NewH - (HtoR + NewD);
    R = R + NewR;
    D = D + NewD;
    
    t=t+1;
end

TotalDeaths = SumD(t-1)
DaysLockdown = sum(lockdown)
DaysCOVID = t

%plotted on a logarithmic y scale
figure;

plot(SumS, 'Color', '#377eb8', 'LineWidth',1.5, 'DisplayName','Susceptible'); hold on
plot(SumE, 'Color', '#62466B', 'LineWidth',1.5, 'DisplayName','Exposed');
plot(SumI, 'Color', '#e41a1c', 'LineWidth',1.5, 'DisplayName','Infected');
plot(SumA, 'Color', '#ff7f00', 'LineWidth',1.5, 'DisplayName','Asymptomatic');
plot(SumH, 'Color', '#984ea3', 'LineWidth',1.5, 'DisplayName','Hospitalised'); %has slight problem
plot(SumR, 'Color', '#4daf4a', 'LineWidth',1.5, 'DisplayName','Recovered');
plot(SumD, 'k','LineWidth',1.5,'DisplayName','Deceased');

figure(2);
plot(lockdown, 'Displayname', 'Lockdown');hold on
plot(mask, 'Displayname', 'Mask')
ylim([0 1.1])
xlim([0 t+5])

% figure(2);
% plot(lockdown, 'Displayname', 'Lockdown');hold on
% plot(mask, 'Displayname', 'Mask')
% ylim([0 1.1])
% xlim([0 t+5])

grid on
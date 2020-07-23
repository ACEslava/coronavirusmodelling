clear all
close all
clc
% Parameters
%load FB2404
%A=ones(2404)-eye(2404);

%Creates the initial graph
Asiz = 5000;
movement = 1; %Initial movement is 1*normal movement
positionx = unifrnd(-1,1,Asiz,1);
positiony = unifrnd(-1,1,Asiz,1);
node1=positionx;
node2=positiony ;
node1(node1>1) = -1;
node1(node1<-1) = 1;
node2(node2>1) = -1;
node2(node2<-1) = 1;

proximity = 0.1;
[A,Distance] = deal(zeros(Asiz,Asiz));
for j=1:Asiz
    for i=1:Asiz
        Distance(i,j) = sqrt((node1(i)-node1(j))^2+(node2(i)-node2(j))^2);
        A(i,j) = Distance(i,j) < proximity;
    end
end

%Parameters
alpha = 0.09; 
beta = 0.1;
gamma = 0.08; 
delta = 0.7; 
theta = 0.99999820042;
h = 0.2033;
omega = 0.01;
Seeds = 1;

%Initialise variables
E = zeros(Asiz,1); %assigns initial seeds
for seed = 1:Seeds
    E(ceil(Asiz*rand),1) = 1;
end

S = ones(Asiz,1) - E; %assigns susceptible population
[I,As,R,H,D] = deal(zeros(Asiz,1)); %all other populations are 0 (no vaccines)
[NewS,NewE,NewI,NewA,NewR,NewH,NewD] = deal(zeros(Asiz,1));
[SumS,SumE,SumI,SumA,SumR,SumH,SumD] = deal(zeros(1000,1));


t = 1;
ACurrent = A;
[mask,lockdown] = deal(zeros(200,1));
[lockdowndays,zerodays] = deal(0);
[FullLockdown,RestrictedLockdown,NoLockdown] = deal(0);
HRandom = rand(Asiz,1);

%Iterations of infection
while (sum(E) + sum(As) + sum(I) > 0) 

    %Randomly moves nodes based on movement multiplier (0.6x movement, 2x movement, etc)
    node1=positionx+movement*(unifrnd(-0.1,0.1,Asiz,1));
    node2=positiony+movement*(unifrnd(-0.1,0.1,Asiz,1));
    
    node1(node1>1) = -1;
    node1(node1<-1) = 1;
    node2(node2>1) = -1;
    node2(node2<-1) = 1;
    

    for j=1:Asiz
        for i=1:Asiz
            Distance(i,j) = sqrt((node1(i)-node1(j))^2+(node2(i)-node2(j))^2);
            A(i,j) = Distance(i,j) < proximity;
        end
    end
    positionx = node1;
    positiony = node2;
    
    %Mask wearers
    Masks = zeros(Asiz,1);
    for maskcount=1:floor(Asiz*0.5) %mask wearing percentage
        Masks(ceil(Asiz*rand)) = 1; %assigns random people to be mask wearers according to percentage defined in for loop
    end

    %Population headcounts
    SumS(t) = sum(S);
    SumE(t) = sum(E);
    SumA(t) = sum(As);
    SumI(t) = sum(I);
    SumH(t) = sum(H);
    SumR(t) = sum(R);
    SumD(t) = sum(D);

    if SumI(t) == 0
        zerodays = zerodays+1;
    else
        zerodays = 0;
    end
    
    %Lockdown Strategy
    if (SumI(t)/Asiz) > 0.005 %Full Lockdown
        proximity = 0.03;
        movement = 0.4;
        lockdown(t) = 1;
    elseif zerodays > 14 %None
        proximity = 0.1;
        movement = 1;
        lockdown(t) = 0;
        lockdowndays = 0;
    elseif t == 1
        lockdown (1) = 0;
    else
        lockdown(t) = lockdown(t-1);
    end
    
    if lockdowndays > 30 %Restricted Lockdown
        proximity = 0.05; 
        movement = 0.8;
        lockdown(t) = 0.5;
        lockdowndays = 0;
    end
    
    if lockdown(t) == 1
        FullLockdown = FullLockdown + 1;
        lockdowndays = lockdowndays + 1;
    elseif lockdown(t) == 0.5
        RestrictedLockdown = RestrictedLockdown + 1;
    else
        NoLockdown = NoLockdown + 1;
    end
    
    %Mask Strategy
    if (SumI(t)/Asiz) > 0.00000094285
        mask(t) = 1;
    elseif zerodays > 21
        mask(t) = 0;
    elseif t==1
        mask(1) =0;
    else
        mask(t) = mask(t-1);
    end

    %Mask wearing neighbors
    INeighbors = A*(or(I,As)); %vector of infected neighbors
    INeighbors_Mask = A*(or(I,As)).*Masks;
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
    %NewH = and((IRandom<0.1), boolean(I));
    age85 = and (boolean(I), and ((HRandom < 0.02), (IRandom < 0.03)));
    age75 = and (boolean(I), and (and ((HRandom > 0.02), (HRandom < 0.06)), (IRandom < 0.02)));
    age65 = and (boolean(I), and(and ((HRandom > 0.06), (HRandom < 0.16)), (IRandom < 0.01)));
    ageother = and (boolean(I), and ((HRandom > 0.16), (IRandom < 0.001)));
    NewH = or (age85, age75);
    NewH = or (NewH, age65);
    NewH = or (NewH, ageother);


    %I or H or A to R
    %R had to differently programmed b/c it has multiple in edges
    ItoR = and((IRandom>1-delta), boolean(I));
    HtoR = and((HRandom<theta), boolean(H));
    AtoR = and((rand(Asiz,1)<delta), boolean(As));
    NewR = or(or(ItoR, HtoR), AtoR);
    
    %H to D
    anotherrand = rand(Asiz, 1);
    age85 = and ((anotherrand < 0.03), and ((HRandom < 0.02), boolean(H)));
    age75 = and ((anotherrand < 0.02), and (and ((HRandom > 0.02), (HRandom < 0.06)), boolean(H)));
    age65 = and ((anotherrand < 0.01), and (and ((HRandom > 0.06), (HRandom < 0.16)), boolean(H)));
    ageother = and ((anotherrand < 0.001), and ((HRandom > 0.16), boolean(H)));
    
    NewD = or (age85, age75);
    NewD = or (NewD, age65);
    NewD = or (NewD, ageother);
    NewD = and (not(NewR), NewD);

    %Update indicators
    S = S - NewE;
    E = E + NewE - (NewI + NewA);
    I = I + NewI - ItoR - NewH;
    As = As + NewA - AtoR;
    H = H + NewH - (HtoR + NewD);
    R = R + NewR;
    D = D + NewD;
    t
    t=t+1;
end
SumS(t:1000) = [];
SumE(t:1000) = [];
SumI(t:1000) = [];
SumR(t:1000) = [];
SumA(t:1000) = [];
SumH(t:1000) = [];
SumD(t:1000) = [];

TotalDeaths = SumD(t-1)
DaysLockdown = sum(lockdown)
DaysCOVID = t
FullLockdown
RestrictedLockdown
NoLockdown

%plotted on a logarithmic y scale
figure;
semilogy(SumS, 'Color', '#377eb8', 'LineWidth',1.5, 'DisplayName','Susceptible'); hold on
semilogy(SumE, 'Color', '#62466B', 'LineWidth',1.5, 'DisplayName','Exposed');
semilogy(SumI, 'Color', '#e41a1c', 'LineWidth',1.5, 'DisplayName','Infected');
semilogy(SumA, 'Color', '#ff7f00', 'LineWidth',1.5, 'DisplayName','Asymptomatic');
semilogy(SumH, 'Color', '#984ea3', 'LineWidth',1.5, 'DisplayName','Hospitalised');
semilogy(SumR, 'Color', '#4daf4a', 'LineWidth',1.5, 'DisplayName','Recovered');
semilogy(SumD, 'k','LineWidth',1.5,'DisplayName','Deceased');

figure(2);
plot(lockdown, 'Displayname', 'Lockdown')
ylim([0 1.1])
xlim([0 t+5])

figure(3);
plot(mask, 'Displayname', 'Mask')
ylim([0 1.1])
xlim([0 t+5])

% figure(2);
% plot(lockdown, 'Displayname', 'Lockdown');hold on
% plot(mask, 'Displayname', 'Mask')
% ylim([0 1.1])
% xlim([0 t+5])

grid on

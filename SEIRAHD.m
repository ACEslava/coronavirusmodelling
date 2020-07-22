clear all
close all
clc
% Parameters
%load FB2404
%A=ones(2404)-eye(2404);

%Creates the initial graph
Asiz = 2404;
movement = 1; %Initial movement is 1*normal movement
positionx = unifrnd(-1,1,Asiz,1);
positiony = unifrnd(-1,1,Asiz,1);
node1=positionx;
node2=positiony ;
node1(node1>1) = -1;
node1(node1<-1) = 1;
node2(node2>1) = -1;
node2(node2<-1) = 1;

proximity = 0.05;
A = zeros(Asiz,Asiz);
for i=1:Asiz
    for j=1:Asiz
        Distance(i,j) = sqrt((node1(i)-node1(j))^2+(node2(i)-node2(j))^2);
        A(i,j) = Distance(i,j) < proximity;
    end
end

%Parameters
alpha = 0.09; %Given
beta = 0.3; %Taken from SEIRDLockdown_NewVersion.m
gamma = 0.08; %(alpha/0.66) = (gamma/0.33)
delta = 0.982; %Given
theta = 0.99999820042; %Given
h = 0.2033; %(delta/0.9) = (h/0.1)
omega = 0.01; %(theta/0.9) = (omega/0.1)
Seeds = 100;

%Initialise variables
E = zeros(Asiz,1); %assigns initial seeds
for seed = 1:Seeds
    E(ceil(Asiz*rand),1) = 1;
end

S = ones(Asiz,1) - E; %assigns susceptible population
[I,As,R,H,D] = deal(zeros(Asiz,1)); %all other populations are 0 (no vaccines)
[NewS,NewE,NewI,NewA,NewR,NewH,NewD] = deal(zeros(Asiz,1));



t = 1;
ACurrent = A;
[mask,lockdown] = deal(zeros(200,1));

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
    

    A = zeros(Asiz,Asiz);
    for i=1:Asiz
        for j=1:Asiz
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

    %Lockdown Strategy, change proximity and movement here. I just picked random values
    if sum(lockdown) > 60
        proximity = 0.03; 
        movement = 0.7;
        lockdown(t) = 0.5;
    elseif SumI(t)/Asiz > 0.02 %Upper threshold for lockdown
        proximity = 0.01;
        movement = 0.3;
        lockdown(t) = 1;
    elseif (SumI(t)/Asiz) < 0.0005 %Lower threshold for lockdown
        proximity = 0.05;
        movement = 1;
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
    age85 = and ((IRandom < 0.03), and ((HRandom < 0.02), boolean(I)));
    age75 = and ((IRandom < 0.02), and(and ((HRandom > 0.02), (HRandom < 0.06)),boolean(I)));
    age65 = and ((IRandom < 0.01), and(and ((HRandom > 0.06), (HRandom < 0.16)),boolean(I)));
    ageother = and ((IRandom < 0.001), and ((HRandom > 0.16), boolean(I)));
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

TotalDeaths = SumD(t-1)
DaysLockdown = sum(lockdown)
DaysCOVID = t

%plotted on a logarithmic y scale
figure;
plot(SumS/Asiz, 'Color', '#377eb8', 'LineWidth',1.5, 'DisplayName','Susceptible'); hold on
plot(SumE/Asiz, 'Color', '#62466B', 'LineWidth',1.5, 'DisplayName','Exposed');
plot(SumI/Asiz, 'Color', '#e41a1c', 'LineWidth',1.5, 'DisplayName','Infected');
plot(SumA/Asiz, 'Color', '#ff7f00', 'LineWidth',1.5, 'DisplayName','Asymptomatic');
plot(SumH/Asiz, 'Color', '#984ea3', 'LineWidth',1.5, 'DisplayName','Hospitalised'); %has slight problem
plot(SumR/Asiz, 'Color', '#4daf4a', 'LineWidth',1.5, 'DisplayName','Recovered');
plot(SumD/Asiz, 'k','LineWidth',1.5,'DisplayName','Deceased');

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

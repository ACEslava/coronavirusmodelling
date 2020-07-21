% SEIRD model

clear all                   % clear the workspace (memory)
close all                   % close all previous plots

% Input parameters
load FB2404                 % loads the adjacency matrix A
N=size(A,1);
%A=ones(2404)-eye(2404);
beta=20/2404;             % spreading rate
alpha=1/5;                  % rate E->I
delta=1/10;                  % rate I->R
omega=(0.10/0.95)*delta;    % rate I->D (assume 10% mortality rate)

% Initialize the state of the population
Seeds=floor(0.005*N);                    % number of exposed seeds
E=zeros(N,1);       % initial infected seed
for seed=1:Seeds
    RandSeed=ceil(N*rand);
    E(RandSeed)=1;
end
S=ones(N,1)-E;     % initial healthy people
I=zeros(N,1); R=zeros(N,1); D=zeros(N,1);

% Run iterations
t=0;
while (sum(E)+sum(I)>0)

    % Update variables
    t=t+1;
    TotalS(t)=sum(S);
    TotalE(t)=sum(E);
    TotalI(t)=sum(I);
    TotalR(t)=sum(R);
    TotalD(t)=sum(D);

    % Transition from S to E
    NI=A*I;  % vector with number of infected neighbors
    NewE=rand(N,1)<1-(1-beta).^(NI);
    NewE=and(NewE,boolean(S));   % boolean vector indicating susceptible people that transfer into Exposed
        
    % Transition from E to I
    CoinE2I=rand(N,1)<alpha;
    NewI=and(CoinE2I,boolean(E));
        
    % Transition from I to R or D
    RandomNumber=rand(N,1);
    CoinI2R=RandomNumber<delta;
    CoinI2D=RandomNumber>1-omega;
    NewR=and(CoinI2R,boolean(I));
    NewD=and(CoinI2D,boolean(I));
    
    % Update indicator vectors
    S=S-NewE;
    E=E+NewE-NewI;
    I=I+NewI-NewR-NewD;
    R=R+NewR;
    D=D+NewD;
end

% Count final results
TotalDeaths=TotalD(t-1)         % total number of deaths
DaysCOVID=t                     % days before eradication

% Plot results
figure;
plot(TotalS/N,'b'); hold on
plot(TotalR/N,'g')
figure;
plot(TotalE/N,'color',[0.9100    0.4100    0.1700]); hold on
plot(TotalI/N,'r')
plot(TotalD/N,'m')
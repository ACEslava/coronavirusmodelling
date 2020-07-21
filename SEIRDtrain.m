% Learning the parameters of the SEIRD simple model
close all
clear all

SEIRD               % generate the true data
%Itrue=TotalI;       % Evolution of infectious
Horizon=t;
TrainHorizon=floor(0.4*Horizon)
Dtrue=TotalD([1:1:TrainHorizon]);       % Evolution of deaths
CumItrue=CumI([1:1:TrainHorizon]);      % Evolution of new infectious

%%
% Initial guess of parameters

b=0.55;     % spreading rate (corrected by 1/N in the eqn below)
a=0.3;      % rate E->I
d=0.4;      % rate I->R
w=0.1*d;   % rate I->D

Seed=10;    % number of initially exposed
N=2404;     % number of individuals

gradsteps=1;
stepsize=1e-5;
Stop=0;

while (Stop==0)
% Compute the model error
[Smodel,Emodel,Imodel,Rmodel,Dmodel,CumImodel]=SEIRDsimpleFunc(TrainHorizon-1,2404,Seeds,b,a,d,w);
ModelError(gradsteps)=sum(((CumItrue-CumImodel)/CumItrue(TrainHorizon)).^2+((Dtrue-Dmodel)/Dtrue(TrainHorizon)).^2)/TrainHorizon;

% Run the model with small increments in each parameter
Delta=0.01;
[Sb,Eb,Ib,Rb,Db,CumIb]=SEIRDsimpleFunc(TrainHorizon-1,2404,Seeds,b+Delta,a,d,w);
[Sa,Ea,Ia,Ra,Da,CumIa]=SEIRDsimpleFunc(TrainHorizon-1,2404,Seeds,b,a+Delta,d,w);
[Sd,Ed,Id,Rd,Dd,CumId]=SEIRDsimpleFunc(TrainHorizon-1,2404,Seeds,b,a,d+Delta,w);
[Sw,Ew,Iw,Rw,Dw,CumIw]=SEIRDsimpleFunc(TrainHorizon-1,2404,Seeds,b,a,d,w+Delta);

% Estimate derivatives as fractions
DerErrorb=((sum(((CumItrue-CumIb)/CumItrue(TrainHorizon)).^2+((Dtrue-Db)/Dtrue(TrainHorizon)).^2)/TrainHorizon)-ModelError(gradsteps))/Delta;
DerErrora=((sum(((CumItrue-CumIa)/CumItrue(TrainHorizon)).^2+((Dtrue-Da)/Dtrue(TrainHorizon)).^2)/TrainHorizon)-ModelError(gradsteps))/Delta;
DerErrord=((sum(((CumItrue-CumId)/CumItrue(TrainHorizon)).^2+((Dtrue-Dd)/Dtrue(TrainHorizon)).^2)/TrainHorizon)-ModelError(gradsteps))/Delta;
DerErrorw=((sum(((CumItrue-CumIw)/CumItrue(TrainHorizon)).^2+((Dtrue-Dw)/Dtrue(TrainHorizon)).^2)/TrainHorizon)-ModelError(gradsteps))/Delta;

% Update the values of the parameters (gradient descent)
b=b-stepsize*DerErrorb;
a=a-stepsize*DerErrora;
d=d-stepsize*DerErrord;
w=w-stepsize*DerErrorw;

if gradsteps>2  % Stopping condition when error increases in a gradient step
    Stop=ModelError(gradsteps)-ModelError(gradsteps-1)>0;
end

gradsteps=gradsteps+1
end

b
a
d
w

[Smodel,Emodel,Imodel,Rmodel,Dmodel,CumImodel]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b,a,d,w);

figure; plot(ModelError); title('Evolution of model error during training')
figure; plot(CumI,'r'); hold on; plot(CumImodel,'r'); title('True infections vs Model infections')
plot(CumItrue); hold on; plot(CumImodel); title('True infections vs Model infections')

figure; plot(TotalD,'m'); hold on; plot(Dmodel,'m'); title('True deaths vs Model deaths')
plot(Dtrue); hold on; plot(Dmodel); title('True deaths vs Model deaths')

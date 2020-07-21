% The model approaches reality as the training increments. Displaying the
% True Infections vs Model Infections with a gradstep less than 1000 showed
% a Model Infection that was a significant underestimation of the True
% Infection. Increasing the gradstep to 2000 showed a Model Infection that
% approached the trend of the True Infection.

% The accuracy of the Model Infection can be evaluated using the Model
% Error function in the script that was graphed as "Evolution of Model
% Error During Training." As the Model Error approached 0, the discrepancy
% between the Model Infection and True Infection decreased. This makes
% sense, as the Model Error function evaluates the difference between the true
% infection number and the model infection number at each gradstep. 

%But as the gradstep increased past a critical point (around 1600), the Model
% Error showed a slight increase. This can be attributed to the Gradient
% Descent algorithm that was used. Thus, it is important to take the Model
% Infection for the gradstep that has the least Model Error to get the most
% accurate model of the True Infection.

% Learning the parameters of the SEIRD simple model
close all
clear all

SEIRD       % generate the true data
Itrue=TotalI;
Dtrue=TotalD;
Horizon=t;

% Initial guess of parameters

b=0.55;     % spreading rate (corrected by 1/N in the eqn below)
a=0.3;      % rate E->I
d=0.4;      % rate I->R
w=0.1*d;   % rate I->D

Seed=10;    % number of initially exposed
N=2404;     % number of individuals

gradsteps=2;
NewModelError = 1;

[Smodel,Emodel,Imodel,Rmodel,Dmodel]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b,a,d,w);
[PreviousModelError, ModelError(1)]=deal(sum((Itrue-Imodel).^2)/Horizon);

while (PreviousModelError >= NewModelError) 
    % Compute the model error
    PreviousModelError = ModelError(gradsteps-1);
    [Smodel,Emodel,Imodel,Rmodel,Dmodel]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b,a,d,w);
    [NewModelError, ModelError(gradsteps)]=deal(sum((Itrue-Imodel).^2)/Horizon);
    
    ModelErrorDeaths(gradsteps) = sum((Dtrue - Dmodel).^2)/Horizon;

    % Run the model with small increments in each parameter
    Delta=0.001;
    [Sdummy,Edummy,Ib,Rdummy,Db]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b+Delta,a,d,w);
    [Sdummy,Edummy,Ia,Rdummy,Da]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b,a+Delta,d,w);
    [Sdummy,Edummy,Id,Rdummy,Dd]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b,a,d+Delta,w);
    [Sdummy,Edummy,Iw,Rdummy,Dw]=SEIRDsimpleFunc(Horizon-1,2404,Seeds,b,a,d,w+Delta);

    % Estimate derivatives as fractions
    DerErrorb=((sum((Itrue-Ib).^2)/Horizon)-ModelError(gradsteps))/Delta + ((sum((Dtrue-Db).^2)/Horizon)-ModelErrorDeaths(gradsteps))/Delta;
    DerErrora=((sum((Itrue-Ia).^2)/Horizon)-ModelError(gradsteps))/Delta  + ((sum((Dtrue-Da).^2)/Horizon)-ModelErrorDeaths(gradsteps))/Delta;
    DerErrord=((sum((Itrue-Id).^2)/Horizon)-ModelError(gradsteps))/Delta  + ((sum((Dtrue-Dd).^2)/Horizon)-ModelErrorDeaths(gradsteps))/Delta;
    DerErrorw=((sum((Itrue-Iw).^2)/Horizon)-ModelError(gradsteps))/Delta  + ((sum((Dtrue-Dw).^2)/Horizon)-ModelErrorDeaths(gradsteps))/Delta;

    % Update the values of the parameters (gradient descent)
    stepsize=1e-8;
    b=b-stepsize*DerErrorb;
    a=a-stepsize*DerErrora;
    d=d-stepsize*DerErrord;
    w=w-stepsize*DerErrorw;

    gradsteps=gradsteps+1
end
b
a
d
w
figure; plot(ModelError); title('Evolution of model error during training')
figure; plot(Itrue,'r'); hold on; plot(Imodel,'g'); title('True infections vs Model infections')
figure; plot(Dtrue,'r'); hold on; plot(Dmodel,'g'); title('True deaths vs Model deaths')
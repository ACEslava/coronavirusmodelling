clear all
close all
SEIRAHD
Horizon = t-1;
TrainHorizon = floor(0.3*Horizon);
TrueI=SumI([1:1:TrainHorizon],1).';
TrueR=SumR([1:1:TrainHorizon],1).';

%Initial Parameter guesses
alpha = 0.3; %Given
beta = 0.55; %Taken from SEIRDLockdown_NewVersion.m
gamma = 0.34; %(alpha/0.66) = (gamma/0.33)
delta = 0.4; %Given
theta = 0.3; %Given
h = 0.015; %(delta/0.9) = (h/0.1)
omega = 0.05; %(theta/0.9) = (omega/0.1)
Seeds = 1;

gradsteps = 1;
stepsize = 1e-7;
Stop = 0;
Increment = 0.01;

while (Stop == 0)
    
    %Compute Model Error
    [ModelS,ModelE,ModelI,ModelR,ModelA,ModelH,ModelD] = SEIRAHDsimpleFunc(alpha,beta,gamma,delta,theta,h,omega,Asiz,Seeds,TrainHorizon) ;
    ModelError(gradsteps) = ModelErrorFunc(TrueI, ModelI, TrueR, ModelR)/TrainHorizon;
    
    %Increment Parameters
     [~,~,Ia,Ra,~,~,~] = SEIRAHDsimpleFunc(alpha+Increment,beta,gamma,delta,theta,h,omega,Asiz,Seeds,TrainHorizon);
     [~,~,Ib,Rb,~,~,~] = SEIRAHDsimpleFunc(alpha,beta+Increment,gamma,delta,theta,h,omega,Asiz,Seeds,TrainHorizon);
     [~,~,Ig,Rg,~,~,~] = SEIRAHDsimpleFunc(alpha,beta,gamma+Increment,delta,theta,h,omega,Asiz,Seeds,TrainHorizon);
     [~,~,Id,Rd,~,~,~] = SEIRAHDsimpleFunc(alpha,beta,gamma,delta+Increment,theta,h,omega,Asiz,Seeds,TrainHorizon);
     [~,~,It,Rt,~,~,~] = SEIRAHDsimpleFunc(alpha,beta,gamma,delta,theta+Increment,h,omega,Asiz,Seeds,TrainHorizon);
     [~,~,Ih,Rh,~,~,~] = SEIRAHDsimpleFunc(alpha,beta,gamma,delta,theta,h+Increment,omega,Asiz,Seeds,TrainHorizon);
     [~,~,Iw,Rw,~,~,~] = SEIRAHDsimpleFunc(alpha,beta,gamma,delta,theta,h,omega+Increment,Asiz,Seeds,TrainHorizon);

     %Estimate derivatives
     DErrorA = (ModelErrorFunc(TrueI, Ia, TrueR, Ra)/TrainHorizon - ModelError(gradsteps))/Increment;
     DErrorB = (ModelErrorFunc(TrueI, Ib, TrueR, Rb)/TrainHorizon - ModelError(gradsteps))/Increment;
     DErrorG = (ModelErrorFunc(TrueI, Ig, TrueR, Rg)/TrainHorizon - ModelError(gradsteps))/Increment;
     DErrorD = (ModelErrorFunc(TrueI, Id, TrueR, Rd)/TrainHorizon - ModelError(gradsteps))/Increment;
     DErrorT = (ModelErrorFunc(TrueI, It, TrueR, Rt)/TrainHorizon - ModelError(gradsteps))/Increment;
     DErrorH = (ModelErrorFunc(TrueI, Ih, TrueR, Rh)/TrainHorizon - ModelError(gradsteps))/Increment;
     DErrorW = (ModelErrorFunc(TrueI, Iw, TrueR, Rw)/TrainHorizon - ModelError(gradsteps))/Increment;
     
     %Gradient Descent
     alpha = alpha - stepsize*DErrorA;
     beta = beta - stepsize*DErrorB;
     gamma = gamma - stepsize*DErrorG;
     delta = delta - stepsize*DErrorD;
     theta = theta - stepsize*DErrorT;
     h = h - stepsize*DErrorH;
     omega = omega - stepsize*DErrorW;
     
     %Stopping Condition
     if gradsteps > 2
         Stop=ModelError(gradsteps) - ModelError(gradsteps-1) > 0;
     end
     
     Loss = ModelError(gradsteps)
     gradsteps = gradsteps + 1;
     
end

alpha
beta
gamma
delta
theta
h
omega

[ModelS,ModelE,ModelI,ModelR,ModelA,ModelH,ModelD] = SEIRAHDsimpleFunc(alpha,beta,gamma,delta,theta,h,omega,Asiz,Seeds,Horizon);
figure(5);
plot(ModelS, 'Color', '#377eb8', 'LineWidth',1.5, 'DisplayName','Susceptible'); hold on
plot(ModelE, 'Color', '#62466B', 'LineWidth',1.5, 'DisplayName','Exposed');
plot(ModelI, 'Color', '#e41a1c', 'LineWidth',1.5, 'DisplayName','Infected');
plot(ModelA, 'Color', '#ff7f00', 'LineWidth',1.5, 'DisplayName','Asymptomatic');
plot(ModelH, 'Color', '#984ea3', 'LineWidth',1.5, 'DisplayName','Hospitalised');
plot(ModelR, 'Color', '#4daf4a', 'LineWidth',1.5, 'DisplayName','Recovered');
plot(ModelD, 'k','LineWidth',1.5,'DisplayName','Deceased');


figure(6)
plot(ModelError)

figure(7)
plot(ModelR); hold on;
plot(SumR);

function Error = ModelErrorFunc(True1, Model1, True2, Model2)
     Error = sum((True1 - Model1).^2 + (True2 - Model2).^2);
end
     
     
     
     
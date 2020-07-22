% Asiz=2404;
% A = zeros(Asiz,Asiz);
% for i=1:Asiz
% 	for j=1:Asiz
% 		Distance(i,j) = sqrt((node1(i)-node1(j))^2+(node2(i)-node2(j))^2);
% 		A(i,j) = Distance(i,j) < rad;
% % 		if A(i,j) ==1
% % 			semilogy([node1(i) node1(j)],[node2(i) node2(j)]) 
% 	end
% end

load FB2404
Asiz = size(A,1);

alpha = 0.09; %Given
beta = 0.85; %Taken from SEIRDLockdown_NewVersion.m
gamma = 0.08; %(alpha/0.66) = (gamma/0.33)
delta = 0.982; %Given
theta = 0.99999820042; %Given
h = 0.2033; %(delta/0.9) = (h/0.1)
omega = 0.00000179957; %(theta/0.9) = (omega/0.1)
Seeds = 100;

[SimpleS,SimpleE,SimpleI,SimpleR,SimpleA,SimpleH,SimpleD]=SEIRAHDsimpleFunc(alpha,beta,gamma,delta,theta,h,omega,Asiz,Seeds,200);

semilogy(SimpleS, 'Color', '#377eb8', 'LineWidth',1.5, 'DisplayName','Susceptible'); hold on
semilogy(SimpleE, 'Color', '#62466B', 'LineWidth',1.5, 'DisplayName','Exposed');
semilogy(SimpleI, 'Color', '#e41a1c', 'LineWidth',1.5, 'DisplayName','Infected');
semilogy(SimpleA, 'Color', '#ff7f00', 'LineWidth',1.5, 'DisplayName','Asymptomatic');
semilogy(SimpleH, 'Color', '#984ea3', 'LineWidth',1.5, 'DisplayName','Hospitalised');
semilogy(SimpleR, 'Color', '#4daf4a', 'LineWidth',1.5, 'DisplayName','Recovered');
semilogy(SimpleD, 'k','LineWidth',1.5,'DisplayName','Deceased');
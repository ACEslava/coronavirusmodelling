function [SimpleS,SimpleE,SimpleI,SimpleR,SimpleA,SimpleH,SimpleD]=SEIRAHDsimpleFunc(a,b,g,d,th,h,w,Asiz,Seed,Horizon)


    % Initialize the state of the population
    [SimpleS,SimpleE,SimpleI,SimpleR,SimpleA,SimpleH,SimpleD] = deal(zeros(1,Horizon-1));
    SimpleS(1)=Asiz-Seed; SimpleE(1)=Seed;
    [SimpleI(1),SimpleR(1),SimpleA(1),SimpleH(1),SimpleD(1)] = deal(0);

    % Run iterations
    for t=1:Horizon-1
        SimpleS(t+1)= SimpleS(t) - ((b/Asiz)*SimpleS(t)*(SimpleI(t)+SimpleA(t)));
        SimpleE(t+1)= SimpleE(t) + ((b/Asiz)*SimpleS(t)*(SimpleI(t)+SimpleA(t))) - ((a+g)*SimpleE(t));
        SimpleI(t+1)=SimpleI(t) + a*SimpleE(t) - (h+d)*SimpleI(t);
        SimpleR(t+1)=SimpleR(t) + d*(SimpleA(t)+SimpleI(t)) + th*SimpleH(t);
        SimpleA(t+1)=SimpleA(t) + g*SimpleE(t) - d*SimpleA(t);
        SimpleH(t+1)=SimpleH(t) + h*SimpleI(t) - (w+th)*SimpleH(t);
        SimpleD(t+1)=SimpleD(t) + w*SimpleH(t);
    end
%     figure(2);
%     plot(SimpleS, 'Color', '#377eb8', 'LineWidth',1.5, 'DisplayName','Susceptible'); hold on
%     plot(SimpleE, 'Color', '#62466B', 'LineWidth',1.5, 'DisplayName','Exposed');
%     plot(SimpleI, 'Color', '#e41a1c', 'LineWidth',1.5, 'DisplayName','Infected');
%     plot(SimpleA, 'Color', '#ff7f00', 'LineWidth',1.5, 'DisplayName','Asymptomatic');
%     plot(SimpleH, 'Color', '#984ea3', 'LineWidth',1.5, 'DisplayName','Hospitalised');
%     plot(SimpleR, 'Color', '#4daf4a', 'LineWidth',1.5, 'DisplayName','Recovered');
%     plot(SimpleD, 'k','LineWidth',1.5,'DisplayName','Deceased');
end


close all
clear all
Asiz = 2404;
positionx = unifrnd(-1,1,Asiz,1);
positiony = unifrnd(-1,1,Asiz,1);
for t=1:5
    node1=positionx+(unifrnd(-0.1,0.1,Asiz,1));
    node2=positiony+(unifrnd(-0.1,0.1,Asiz,1));
    
    for k=1:Asiz
        if rand < 0.01
        node1(k) = unifrnd(-0.1,0.1);
        node2(k) = unifrnd(-0.1,0.1);
        end
    end
        
    node1(node1>1) = -1;
    node1(node1<-1) = 1;
    node2(node2>1) = -1;
    node2(node2<-1) = 1;
    figure; plot(node1,node2, 'o'); hold on
%     rad = 0.05;
%     A = zeros(Asiz,Asiz);
%     for i=1:Asiz
%         for j=1:Asiz
%             Distance(i,j) = sqrt((node1(i)-node1(j))^2+(node2(i)-node2(j))^2);
%             A(i,j) = Distance(i,j) < rad;
%         end
%     end
%     positionx = node1;
%     positiony = node2;
    t
end
close all
clear all
Asiz = 2404;
node1=rand(Asiz,1);
node2=rand(Asiz,1);

%figure; plot(node1,node2,'o'); hold on
rad=0.03;

A = zeros(Asiz,Asiz);
for i=1:Asiz
	for j=1:Asiz
		Distance(i,j) = sqrt((node1(i)-node1(j))^2+(node2(i)-node2(j))^2);
		A(i,j) = Distance(i,j) < rad;
% 		if A(i,j) ==1
% 			plot([node1(i) node1(j)],[node2(i) node2(j)]) 
	end
end
spy(A)

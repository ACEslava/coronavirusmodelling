function [S,E,I,R,D]=SEIRDsimpleFunc(Horizon,N,Seed,b,a,d,w)

% Initialize the state of the population

S(1)=N-Seed; E(1)=Seed;
I(1)=0;  R(1)=0; D(1)=0;

% Run iterations
for t=1:Horizon
    S(t+1)=S(t)-(b/N)*S(t)*I(t);
    E(t+1)=E(t)+(b/N)*S(t)*I(t)-a*E(t);
    I(t+1)=I(t)+a*E(t)-(d+w)*I(t);
    R(t+1)=R(t)+d*I(t);
    D(t+1)=D(t)+w*I(t);
end

end

% LA1D_m.m
clc
clear all

% constants
N=20;x0=0.;xN=1.;deltaT=0.01;c=1.;T=0.1;
% object LA1D initialization
LA1D=LinearAdvection1D(N,x0,xN,deltaT,c,T);

x=linspace(LA1D.x0,LA1D.xN,LA1D.N);
% initial value
u0=zeros(1,LA1D.N);
for i=5:10
   u0(i)=1;
end
% plot of initial value
plot(x,u0)
xlabel("x")
ylabel("u")

% calculating solution if CFL<=1
if (LA1D.checkCFL()== true)
	disp(strcat("CFL number is: ",num2str(LA1D.CFL())))
	LA1D.upwindMatrixAssembly()
	for t=0:uint8(LA1D.T/LA1D.deltaT)
		u=LA1D.upwindSolve(u0);
		u0=u;
	end
else
	disp(strcat("CFL number is greater than 1. CFL: ",num2str(LA1D.CFL())))
end

% ploting the last solution
hold on
plot(x,u)
legend({"Initial value",strcat("Solution at t=",num2str(LA1D.T))},'Location','east')
axis([0 1.1  0 1.1]) 
grid on


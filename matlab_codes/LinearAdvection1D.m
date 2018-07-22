	classdef LinearAdvection1D < handle    
		properties
			%  constants
			A;N;x0;xN;deltaT;c;T;
		end
		
		methods
			% Initialization of constants
			function self = LinearAdvection1D(N,x0,xN,deltaT,c,T) 
				self.N = N; self.x0 = x0; self.xN = xN; 
				self.deltaT = deltaT;
				self.c = c; self.T = T;self.A = zeros(N);
			end
			% CFL number funct.
			function outputArg = CFL(self) 
				deltaX= (self.xN - self.x0)/self.N ;
			   outputArg = abs(self.c*self.deltaT/deltaX);
			end
			% check CFL number <=1 or not.
			function flag=checkCFL(self)
			   if (self.CFL()<=1)
					flag=true;
			   else
					flag=false;
			   end 
			end
			% Matrix assembly of LA1D
			function upwindMatrixAssembly(self)
			   alpha_min=min(self.c,0);
			   alpha_max=max(self.c,0);
			   a1=alpha_max*ones(1,self.N - 1);
			   a2=(1+alpha_min-alpha_max)*ones(1,self.N);
			   a3=-alpha_min*ones(1,self.N-1);
			   self.A=diag(a1, -1)+diag(a2, 0)+diag(a3, 1);
			   self.A(1,end)=alpha_max;
			   self.A(end,1)=-alpha_min;
			end
			% Upwind solve
			function x = upwindSolve(self,u0)
				x=transpose(self.A\transpose(u0));
			end
		end
	end

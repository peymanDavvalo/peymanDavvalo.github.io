 
type LinearAdvection1d
    # dataIn  [c,x0,xN,deltaT,T]
    dataIn::Array{Float64,1}     
    N::Int    
    
    Initialize::Function
    CFL::Function
    checkCFL::Function
    upwindMatrixAssembly::Function
    upwindSolve::Function
    
    function LinearAdvection1d()
        self = new()
        
        self.Initialize = function (dataIn::Array{Float64,1}, N::Int)
            self.dataIn = dataIn
            self.N = N
        end
        self.CFL = function () 
            deltaX= (self.dataIn[3] - self.dataIn[2])/self.N
            return abs(self.dataIn[1]*self.dataIn[4]/deltaX)
        end
        self.checkCFL = function ()
            return  self.CFL()<=1 ? true : false
        end

        self.upwindMatrixAssembly = function() 
            alpha_min=min(self.dataIn[1],0)
            alpha_max=max(self.dataIn[1],0)
            a1=[alpha_max for n in 1:self.N-1]
            a2=[1+alpha_min-alpha_max for n in 1:self.N]
            a3=[-alpha_min for n in 1:self.N-1]
            A=Tridiagonal(a1,a2,a3)+zeros(self.N,self.N)
            A[1,end]=alpha_max
            A[end,1]=-alpha_min
            return A
        end
        
        self.upwindSolve = function(u0::Array{Float64,1})
             return self.upwindMatrixAssembly()\u0
        end        
    end 
end

#####################
# Start of the code
#####################

N,x0,xN,deltaT,c,T=20,0.,1.,0.01,1.,0.1
LA1D=LinearAdvection1d()
LA1D.self.Initialize([c;x0;xN;deltaT;T],N)
LA1D.self.checkCFL()

u0 = [0. for n in 1:N];
u0[5:10]=1;

PyPlot.plot(x,u0,label="Initial value")


x=linspace(x0,xN,N);
     
        
if LA1D.self.checkCFL()
    println("CFL number is: ", LA1D.self.CFL())
    for t=0 : floor(Integer,T/deltaT)
        u=LA1D.self.upwindSolve(u0)
        u0=u
    end
else
    println("CFL number is greater than 1. CFL: ", LA1D.self.CFL())    
end

PyPlot.plot(x,u0,label=string("Solution at t=",T))
PyPlot.legend()
PyPlot.grid(linestyle="dotted")


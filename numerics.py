from sympy import *
import numpy as np
from math import *
import matplotlib.pyplot as plt 
class ODE:
    def __init__(self,func):
        self.order=len(func) #order of the differential equation
        self.func=func #list of functions into which the differential equation has been decomposed

    def euler(self,init,x0,high,step):
        #init: list of initial values passed in the same order as functions
        #x0: initial value of the independent variable
        #high: 
        #step: the step size
        x=np.arange(x0,high,step)
        temp=init
        y=np.zeros_like(x)
        for i in range(len(x)):
          j=0
          init=temp
          for func in self.func:
            temp[j]=init[j]+step*func(x[i],*init)
            j+=1 
          y[i]=temp[-1]
        return x,y
    
    def rk4(self,init,x0,high,step):
       #init: list of initial values passed in the same order as functions
        #x0: initial value of the independent variable
        #high: 
        #step: the step size
       x=np.arange(x0,high,step)
       init=np.array(init)
       temp=init.copy()
       y=np.zeros_like(x)
       for i in range(len(x)):
          k=np.zeros_like(init)
          ki=np.zeros_like(init)
          d=np.zeros_like(init)
          poly1=lambda x: (2*(x**3)-9*(x**2)+(13*x))/12
          poly2=lambda x: (-(x**2)+(3*x)+2)/2
          for n in range(4):
            j=0
            for func in self.func:
              k[j]=step*func(x[i]+(step*poly1(n)),*list(temp+(ki*poly1(n))))
              j+=1
            ki=k
            d+=poly2(n)*k
          temp+=(d/6)
          y[i]=temp[-1]
       return x,y
    
class RootFinder:
    def __init__(self,f):
        x,y,z=symbols('x y z')
        self.deriv=lambdify(x,diff(sympify(f),x),"math") #passing the derivative of the given expression
        self.expres=lambdify(x,sympify(f),"math") #passing the expression

    def Bisection(self,a,b,e=0.01):
        #a: the lower bound of the interval
        #b: the upper bound of the interval
        #e: the accuracy
        x,y,z=symbols('x y z')
        if (self.expres(a)<0 and self.expres(b)>0 ):
            low=a
            high=b
        elif self.expres(a)*self.expres(b)>0:
            raise Exception("no root in the given range (function doesn't change sign)")
        else:
            low=b
            high=a
        n=ceil((log(b-a)-log(e))/log(2))
        for i in range(n):
            root=(high+low)/2
            if self.expres(root)>0:
                high=root
            else:
                low=root
        return root
    
    def NewtonRaphson(self,x0):
        #x0: initial guess
        root=x0
        n=1
        while abs(self.expres(root))>0.01:
            if self.deriv(root)==0:
                raise ZeroDivisionError
            root=root-(self.expres(root)/self.deriv(root))
            n+=1
            if n>1000:
                raise Exception("initial guess inaccurate, method did not converge")
        return root
    
class Integral:
    def __init__(self,f):
        self.func=f #passing the function of which the integral has to be found

    def trapezoidal(self,xi,xf,step):
        #x1,xf: the lower and upper bounds of the interval in which the integral has to be found
        #step:step size
        integral=0
        interval=np.arange(xi,xf+step,step)
        for i in range(len(interval)-1):
          integral+=(self.func(interval[i])+self.func(interval[i+1]))
        return (step/2)*integral
    
    def booles(self,xi,xf,step):
        #x1,xf: the lower and upper bounds of the interval in which the integral has to be found
        #step:step size
        if ((xf-xi)/step)%4!=0:
            raise Exception("number of intervals must be a multiple of four")
        x=np.arange(xi,xf+step,step)
        integral=0
        for i in range(0,len(x)-4,4):
          integral+=((7*self.func(x[i]))+(32*self.func(x[i+1]))+(12*self.func(x[i+2]))+(32*self.func(x[i+3]))+(7*self.func(x[i+4])))
        return (2*step/45)*integral  
          


          

            

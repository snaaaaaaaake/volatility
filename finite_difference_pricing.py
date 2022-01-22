import matplotlib.pyplot as plt
from scipy.stats     import norm
from math  			 import *
from numpy 			 import *
from volatilities    import *

##### Generic parameters #####
n       = 1000
p       = 300
K       = 170.0
T       = 0.235455
s_infty = 340.0
dt      = T/n
ds      = s_infty/p
I       = list(range(0,n))
J       = list(range(0,p+1))
I.reverse()

nT = n
pK = p
dT = dt
dK = ds
IT = list(range(0,nT))
JK = list(range(0,pK+1))
IT.reverse()


##### Black-Scholes pricing #####
# Normal cumulative distribution function
def N(x):
	return norm.cdf(x)

def σ_bs(t, s):
	return 0.4	

def V_bs(t, s, T = T, K = K, σ = σ_bs, ε = 1):
	dt = T - t
	d  = exp(-r*dt)
	F  = s/d

	if F == 0:
		F = 0.0000001

	if K == 0:
		K = 0.0000001

	if dt == 0:
		dt = 0.0000001

	d1 = (log(F/K) + 0.5*dt*σ(dt,K)**2)/(σ(dt,K)*sqrt(dt)) if F != 0 else -10000000
	d2 = d1 - σ(dt,K)*sqrt(dt)
	return ε*d*(F*N(ε*d1) - K*N(ε*d2))

##### Financial contracts #####
# vanilla call option payoff
def call_option(s):
	return max(s - K,0)

def put_option(s):
	return max(K - s,0)

def Φ_us(t, s, v, payoff, T):
	return max(v, payoff(s)) if t < T else payoff(s)

def Φ_eu(t, s, v, payoff, T):
	return v if t < T else payoff(s)
##### Market data #####
r = 0.06

#def σ_imp(T, K):
#	return volatility[T,K]

def V(t, s, T, K):
	return V_bs(t=t,s=s,T=T,K=K,ε=1, σ=σ_bs)*exp(r*(T-t))

def σ_loc(T, K):
	if K == 0:
		K = 0.0001

	t = 0
	s = 170
	Ckk = (V(t,s,T,K + 2*dK) + V(t,s,T,K) - 2*V(t,s,T,K+dK))/(dK**2)
	CT  = (V(t,s,T+dT,K) - V(t,s,T,K))/dT
	return sqrt(2*CT/Ckk)/K

def σ(t, s):
	return σ_loc(t,s)*s


##### PDE functions #####
# a, b, c and d are chosen to refer to vt + 0.5*σ*σ.vss + r.s.vs -rv = 0
# corresponding to the SDE dS(t) = r.S(t).dt + σ(t,S(t)).dW(t)
def a(t, s):
	return 1.

def b(t, s):
	return 0.5*σ(t,s)**2

def c(t, s):
	return r*s

def d(t, s):
	return -r

def Ψ(i, j, ε):
	return (b(i,j)/ds + 0.5*ε*c(i,j))/ds

def A(i, j):
	if(j == 0):
		return 0.5*Ψ(i,j,0)
	elif(j == p):
		return d(i,j) - a(i,j)/dt + 0.5*Ψ(i,j,4)
	else:
		return Ψ(i,j,+1)

def B(i, j):
	if(j == 0):
		return -Ψ(i,j,-2)
	elif(j == p):
		return -Ψ(i,j,+2)
	else:
		return d(i,j) - a(i,j)/dt - 2*Ψ(i,j,0)

def C(i, j):
	if(j == 0):
		return d(i,j) - a(i,j)/dt + 0.5*Ψ(i,j,-4)
	elif(j == p):
		return 0.5*Ψ(i,j,0)
	else:
		return Ψ(i,j,-1)

def λ(i, j):
	return -a(i,j)/dt

# solve a.vt + b.vxx + c.vx + d.v = 0
def solve(payoff = put_option, Φ = Φ_eu, T = T):
	# Resolution of the PDE
	Vi = zeros((p+1,1))
	for i in I:
		λi = diag([λ(i,j) for j in J])
		Mi = zeros((p+1,p+1))

		# Boundary condition
		Vi_ = array([ [Φ((i+1)*dt,j*ds,Vi[j,0],payoff,T)] for j in J ])

		for j in J:
			if(j == 0):
				Mi[0,0] = C(i,0)
				Mi[0,1] = B(i,0)
				Mi[0,2] = A(i,0)
			elif(j == p):
				Mi[p,p-2] = C(i,p)
				Mi[p,p-1] = B(i,p)
				Mi[p,p]   = A(i,p)
			else:
				Mi[j,j-1] = C(i,j)
				Mi[j,j]   = B(i,j)
				Mi[j,j+1] = A(i,j)

		# Back propagation
		Mi_inv = linalg.inv(Mi)
		Vi = Mi_inv.dot(λi)
		Vi = Vi.dot(Vi_) 

	return Vi




##### Main #####
#payoff = call_option
#
#Veu = solve(Φ = Φ_eu, payoff = payoff) # PDE pricing function of a european payoff with constant vol 
#Vus = solve(Φ = Φ_us, payoff = payoff)          # PDE pricing function of an american payoff with constant vol
#Vbs = [V_bs(0,j*ds,ε=1 if payoff==call_option else -1) for j in J]   # Black scholes pricing function with constant vol
##Vbss = [V_bs(0,j*ds,ε=-1, σ=σ_imp) for j in J]   # Black scholes pricing function with constant vol
#
#
#plt.plot([j*ds for j in J],list(Veu[:,0]),'b')
#plt.plot([j*ds for j in J],list(Vus[:,0]),'r')
#plt.plot([j*ds for j in J],Vbs,'y')
#plt.plot([j*ds for j in J],list(array([ [payoff(j*ds)] for j in J ])[:,0]),'g')
#
#plt.axis([0,300,0,200])
#plt.show()
#
#print("Solved !")



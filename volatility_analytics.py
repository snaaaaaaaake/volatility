from finite_difference_pricing import *



def σ_(t, s):
	return 0.4

# undiscounted price

def pd(t, s, T, K, σ = σ_):

	µ = log(s) + (r - 0.5*σ(t,s)**2)*(T-t)
	sg = σ(t,s)*sqrt(T-t)

	if(K==0):
		K = 0.0000001

	A = -(log(K) - µ)**2
	B = 2*sg**2
	C = K*sg*sqrt(2*pi)

	return exp(A/B)/C

(t,s,T) = (0,170.,2.0)

C   = [V(t,s,T,k*dK) for k in JK]
Ck  = [(V(t,s,T,(k+1)*dK) - V(t,s,T,k*dK))/dK for k in JK] 
Ckk = [(V(t,s,T,(k+2)*dK) + V(t,s,T,k*dK) - 2*V(t,s,T,(k+1)*dK))/(dK**2) for k in JK] 
CT  = [(V(t,s,T+dT,k*dK) - V(t,s,T,k*dK))/dT for k in JK]
pp  = [pd(t,s,T,k*dK) for k in JK]

def σ_loc(t, s, T, K):
	if K == 0:
		K = 0.0001

	Ckk = (V(t,s,T,K + 2*dK) + V(t,s,T,K) - 2*V(t,s,T,K+dK))/(dK**2)
	CT  = (V(t,s,T+dT,K) - V(t,s,T,K))/dT
	return sqrt(2*CT/Ckk)/K

#plt.plot(JK,C,'b')
#plt.plot(JK,Ck,'r')
#plt.plot(JK,Ckk,'g')
#plt.plot(JK,pp,'y')
plt.plot([k*dK for k in JK if k > 50],[σ_loc(0,170.,2.0,k*dK) for k in JK if k > 50],'b')
plt.plot([k*dK for k in JK if k > 50],[σ_loc(0,100.,2.0,k*dK) for k in JK if k > 50],'g')

#plt.plot(JK,CT,'r')

plt.show()
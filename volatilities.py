import pandas as pd
from matplotlib import pyplot as plt

#path = "/Users/neo/Ensimag/PhD/aapl_surface.csv"
#surface = pd.read_csv(path, sep=',')
#
#T = 'EXPIRY'
#K = 'STRIKE'
#σ = 'VOLATILITY'
#
#dT = surface[T]
#T = list(dT.unique())
#dK = surface[K]
#K = list(dK.unique())
#dσ = surface[σ]
#
#volatility = {(T,K):σ for (T,K,σ) in zip(list(dT),list(dK),list(dσ))}

def σ_imp(T, K):
#	return volatility[T,K]
	return 0.2


def plot_function(x, function, color = 'b'):
	plt.plot(x, [function(xx) for xx in x], color)
	plt.show()
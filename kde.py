from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt  
import numpy as np   


# ╦╔═┌─┐┬─┐┌┐┌┌─┐┬    ╔╦╗┌─┐┌┐┌┌─┐┬┌┬┐┬ ┬  ╔═╗┌─┐┌┬┐┬┌┬┐┌─┐┌┬┐┬┌─┐┌┐┌
# ╠╩╗├┤ ├┬┘│││├┤ │     ║║├┤ │││└─┐│ │ └┬┘  ║╣ └─┐ │ ││││├─┤ │ ││ ││││
# ╩ ╩└─┘┴└─┘└┘└─┘┴─┘  ═╩╝└─┘┘└┘└─┘┴ ┴  ┴   ╚═╝└─┘ ┴ ┴┴ ┴┴ ┴ ┴ ┴└─┘┘└┘


x = [1,2,-1,-1,-1,4,5]
density = gaussian_kde(x)
xgrid = np.linspace(min(x), max(x), len(x))   
plt.plot(xgrid, density(xgrid))
plt.show()

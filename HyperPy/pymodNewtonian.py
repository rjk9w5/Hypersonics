import math

def newtonian(theta):
    return 2 * math.sin(theta)

def modifiedNewtonian(theta):
    cpmax = 2
    return cpmax * math.sin(theta)

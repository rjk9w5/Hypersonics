import math
import scipy.optimize as sciO

class cpg:

    class obs:

        def __init__(self, delta, M, gamma=1.4):
            self.delta = delta
            self.M = M
            self.gamma = gamma

        def shockAngle(self):

            # Estimate beta for a given mach number and deflection angle
            B0 = (60 / self.M * math.pi / 180 + self.delta)

            self.beta = sciO.newton(self.__dbm__, B0)
            return self.beta


        def __dbm__(self,beta):
            return 2 * (1 / (math.tan(beta))) * ((self.M ** 2 * math.sin(beta) ** 2 - 1) /
                    (self.M ** 2 * (self.gamma + math.cos(2 * beta)) + 2)) - math.tan(self.delta)

    class nsw:

        def __init__(self, M, gamma=1.4):
            self.M = M
            self.gamma = gamma

        def prat(self):
            return

        def Trat(self):
            return

        def rrat(self):
            return

class equillibrium:

    def __init__(self):
        self.empty=1
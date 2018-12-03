
class galaxy:

    def __init__(self, theta=None, phi=None, z=None, childs=None):
        self._theta  = theta
        self._phi    = phi
        self._z      = z
        self._childs = childs

    @property
    def theta(self):
        return self._theta
    @theta.setter
    def theta(self, theta):
        self._theta = theta

    @property
    def ra(self):
        return 90.-self._phi

    @property
    def dec(self):
        return self._theta
    
    @property
    def phi(self):
        return self._phi
    @phi.setter
    def phi(self, phi):
        self._phi = phi

    @property
    def z(self):
        return self._z
    @z.setter
    def z(self, z):
        self._z = z

    @property
    def childs(self):
        return self._childs
    @childs.setter
    def childs(self, thetas, phis, zs):
        num_childs = len(ras)
        for i in range(num_childs):
            self._childs[i] = galaxy(thetas[i], phis[i], zs[i])

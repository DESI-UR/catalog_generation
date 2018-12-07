
TYPE_LUT = {0: "center", 1: "rim", 2: "flat", 3: "clump_center", 4: "clump_flat"}

class galaxy:

    def __init__(self, name=None, theta=None, phi=None, r=None, childs=[], parent=None, TYPE=None):
        self._name   = name
        self._theta  = theta
        self._phi    = phi
        self._r      = r
        self._TYPE   = TYPE
        self._childs = childs
        self._parent = parent

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = name
        
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
    def r(self):
        return self._r
    @r.setter
    def r(self, r):
        self._r = r

    @property
    def TYPE(self):
        return self._TYPE
    @TYPE.setter
    def TYPE(self, TYPE):
        self._TYPE = TYPE
        
    @property
    def childs(self):
        return self._childs
    @childs.setter
    def childs(self, childs):
        self._childs = childs

    @property
    def parent(self):
        return self._parent
    @parent.setter
    def parent(self, parent):
        self._parent = parent
        
    def print_info(self, print_childs=False):
        print("Galaxy ({})".format(TYPE_LUT[self._TYPE]))
        print("  Name:   {}".format(self.name))
        print("  RA:     {}".format(self.ra[0]))
        print("  Dec:    {}".format(self.dec[0]))
        print("  Dist:   {}".format(self.r))
        if self.parent is not None:
            print("  Prnt:   {}".format(self.parent))
        """
        try:
            print("  Childs: {}".format(len(self.childs)))
            if print_childs:
                for i in range(len(self.childs)):
                    self.childs[i].print_info(print_childs=True)
        except:
            return
        """

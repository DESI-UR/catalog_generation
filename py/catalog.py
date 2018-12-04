
class catalog:

    def __init__(self):
        self._centers = None
        self._rims    = None
        self._flats   = None
        self._clumps_center = None
        self._clumps_flat   = None

    @property
    def centers(self):
        return self._centers
    @centers.setter
    def centers(self, centers):
        self._centers = centers

    @property
    def rims(self):
        return self._rims
    @rims.setter
    def rims(self, rims):
        self._rims = rims

    @property
    def flats(self):
        return self._flats
    @flats.setter
    def flats(self, flats):
        self._flats = flats

    @property
    def clumps_center(self):
        return self._clumps_center
    @clumps_center.setter
    def clumps_center(self, clumps_center):
        self._clumps_center = clumps_center

    @property
    def clumps_flat(self):
        return self._clumps_flat
    @clumps_flat.setter
    def clumps_flat(self, clumps_flat):
        self._clumps_flat = clumps_flat

    
        

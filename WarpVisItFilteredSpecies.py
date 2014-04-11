import numpy
import warp

#############################################################################
class WarpVisItFilteredSpecies(object):
    """
    Interface to Warp Species class that supports filtering
    Construct it with a species instance and a filter function
    that generates a set of ids selecting the desired particles
    declared like: def filterFunc(x,y,z,ux,uy,uz): ... return ids
    """
    #------------------------------------------------------------------------
    def __init__(self,species,filterFunc,filterName):
        """ """
        self.__Species = species
        self.__Filter = filterFunc
        self.__Ids = None
        self.__Key = -1
        self.name = species.name + filterName
        self.type = species.type
        return

    #------------------------------------------------------------------------
    def Update(self, key):
        """
        Generate the ids, by calling user provided function. Does
        something only if key has changed since the last invokation.
        """
        if (self.__Ids is None) or (self.__Key != key):
            self.__Ids = self.__Filter(self.__Species.getx(gather=0),
                    self.__Species.gety(gather=0),
                    self.__Species.getz(gather=0),
                    self.__Species.getux(gather=0),
                    self.__Species.getuy(gather=0),
                    self.__Species.getuz(gather=0))
            self.__Key = key
        return

    #------------------------------------------------------------------------
    def Filter(self, data):
        """
        Return filtered data
        """
        self.Update(warp.warp.top.it)
        if (len(self.__Ids) > 0):
            return numpy.take(data,self.__Ids)
        else:
            return numpy.array([])

    #------------------------------------------------------------------------
    def getpid(self,**kwargs): return self.Filter(self.__Species.getpid(**kwargs))
    def getn(self,**kwargs): return self.Filter(self.__Species.getn(**kwargs))
    def getx(self,**kwargs): return self.Filter(self.__Species.getx(**kwargs))
    def gety(self,**kwargs): return self.Filter(self.__Species.gety(**kwargs))
    def getz(self,**kwargs): return self.Filter(self.__Species.getz(**kwargs))
    def getw(self,**kwargs): return self.Filter(self.__Species.getw(**kwargs))
    def getr(self,**kwargs): return self.Filter(self.__Species.getr(**kwargs))
    def gettheta(self,**kwargs): return self.Filter(self.__Species.gettheta(**kwargs))
    def getvtheta(self,**kwargs): return self.Filter(self.__Species.getvtheta(**kwargs))
    def getetheta(self,**kwargs): return self.Filter(self.__Species.getetheta(**kwargs))
    def getbtheta(self,**kwargs): return self.Filter(self.__Species.getbtheta(**kwargs))
    def getvx(self,**kwargs): return self.Filter(self.__Species.getvx(**kwargs))
    def getvy(self,**kwargs): return self.Filter(self.__Species.getvy(**kwargs))
    def getvz(self,**kwargs): return self.Filter(self.__Species.getvz(**kwargs))
    def getvr(self,**kwargs): return self.Filter(self.__Species.getvr(**kwargs))
    def getux(self,**kwargs): return self.Filter(self.__Species.getux(**kwargs))
    def getuy(self,**kwargs): return self.Filter(self.__Species.getuy(**kwargs))
    def getuz(self,**kwargs): return self.Filter(self.__Species.getuz(**kwargs))
    def getex(self,**kwargs): return self.Filter(self.__Species.getex(**kwargs))
    def getey(self,**kwargs): return self.Filter(self.__Species.getey(**kwargs))
    def getez(self,**kwargs): return self.Filter(self.__Species.getez(**kwargs))
    def geter(self,**kwargs): return self.Filter(self.__Species.geter(**kwargs))
    def getbx(self,**kwargs): return self.Filter(self.__Species.getbx(**kwargs))
    def getby(self,**kwargs): return self.Filter(self.__Species.getby(**kwargs))
    def getbz(self,**kwargs): return self.Filter(self.__Species.getbz(**kwargs))
    def getbr(self,**kwargs): return self.Filter(self.__Species.getbr(**kwargs))
    def getxp(self,**kwargs): return self.Filter(self.__Species.getxp(**kwargs))
    def getyp(self,**kwargs): return self.Filter(self.__Species.getyp(**kwargs))
    def getzp(self,**kwargs): return self.Filter(self.__Species.getzp(**kwargs))
    def getrp(self,**kwargs): return self.Filter(self.__Species.getrp(**kwargs))
    def gettp(self,**kwargs): return self.Filter(self.__Species.gettp(**kwargs))
    def getgaminv(self,**kwargs): return self.Filter(self.__Species.getgaminv(**kwargs))
    def getweights(self,**kwargs): return self.Filter(self.__Species.getweights(**kwargs))
    def getke(self,**kwargs): return self.Filter(self.__Species.getke(**kwargs))
    def getrank(self,**kwargs): return self.Filter(self.__Species.getrank(**kwargs))

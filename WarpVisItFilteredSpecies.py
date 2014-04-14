import numpy
import warp
import inspect
import types

#############################################################################
class WarpVisItFilteredSpecies(object):
    """
    Interface to Warp Species class that supports filtering
    Construct it with a species instance and a filter function
    that generates a set of ids selecting the desired particles
    declared like: def filterFunc(x,y,z,ux,uy,uz): ... return ids
    """
    #------------------------------------------------------------------------
    def __init__(self,species,filterFunc,filterName,**kwargs):
        """
        Initalize the WarpVisItFilteredSpecies object:
        
        Parameters: 
        
            species : The species object to which filtering is applied
            filterFunc : The filter function to be applied
            filterName : The name of the filter
            kwargs : Additional parameters to be passed to the filter function
            
        Attributes:
        
            __Species : The species object to which the filter is applied to
            __FilterFunc : The filter funtion
            __FilterArgs : Additional parameters for the filter function
            __IDs : List of selected ids to avoid repeated evalutation of the filter.
            __Key : Timestep key when the filter was last evaluated. Used to check __Ids.
            name : Name for the filtered species, defined by filterName and species.name
            type : Type of the species defined by species.type.
            
        
        """
        # Check if the filterfunction is valid
        if 'species' not in inspect.getargspec(filterFunc).args:
            raise KeyError("filterfunc must have species as input argument.")
        #Set __FilterFunc and __FilterArgs using the approbirate set functions
        self.SetFilter(filterFunc)
        self.SetFilterArgs(**kwargs)
        # Define other object variables
        self.__Species = species
        self.__Ids = None
        self.__Key = -1
        self.name = species.name + filterName
        self.type = species.type
        
        return

    #------------------------------------------------------------------------
    def Update(self, key):                                      
        """
        Generate the ids, by calling the user provided filter function. Ids are
        only updated if the key has changed or Ids are not available.
        
        Parameters:
        
            key : Timestep key used to check whether the existing IDs are valid.
             
        """
        if (self.__Ids is None) or (self.__Key != key):
            self.__Ids = self.__FilterFunc(self.__Species, **self.__FilterArgs)
            self.__Key = key
        return

    
    #------------------------------------------------------------------------
    def SetFilter(self, filterFunc):         
        """
        Change the filter function used.
        
        Parameters:
        
            filterFunc : The new filter function to be used.
        
        """
        # Determine whether we need to bind the filter function to the object
        if inspect.getargspec(filterFunc).args[0] == 'self':
            self.__FilterFunc = types.MethodType(filterFunc, self)
        else:
            self.__FilterFunc = filterFunc
        self.__Ids = None
    
    #------------------------------------------------------------------------
    def SetFilterArgs(self, **kwargs):                         
        """
        Change filter parameters
        
        Parameters:
        
            **kwargs: Define any keyword parameters used by the filter function
        
        """
        self.__FilterArgs = kwargs
        self.__Ids = None

    #------------------------------------------------------------------------
    def Filter(self, data):
        """
        Return filtered data. Filter the given data array. Use the provided 
        get functions to access common variables associated with the species,
        e.g, getx, gety, getz etc. 
        
        Parameters:
        
            data : The data array to be filtered
             
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







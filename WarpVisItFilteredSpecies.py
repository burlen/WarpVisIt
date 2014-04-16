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




#############################################################################
class WarpVisItSpeciesFilters(object):
    """
    Class providing a collection of common filters used with the
    WarpVisItFilteredSpecies class to generate filtered species.
    """
    
    #----------------------------------------------------------------------------
    @classmethod
    def GetAvailableFilters(cls):
        """
        Get dictionary of available filter methods defined by 
        the WarpVisItSpeciesFilters class
        """
        return {'ThresholdFilter': cls.ThresholdFilter,
                'AccumulativeThresholdFilter': cls.AccumulativeThresholdFilter,
                'PidFilter': cls.PidFilter,
                'HaloFilter': cls.HaloFilter}
        

    #----------------------------------------------------------------------------
    @classmethod
    def CreateFilteredSpecies(cls,species,filterType,filterName,**kwargs):
        """
        Convenience factory method to create new filtered species for 
        the given species, using the given filter type.
           
        Paramters:
        
            species : The species object to which filtering is applied
            filterType : The type of filter to be used. One of AvailableFilters
        
        Optional Parameters:
        
            kwargs : Additional parameters specific to the filter functions
        
        """
        filterFunc = cls.GetAvailableFilters()[filterType]
        return WarpVisItFilteredSpecies(species=species,
                                        filterFunc=filterFunc,
                                        filterName=filterName,
                                        **kwargs)
    
    
    #----------------------------------------------------------------------------
    @staticmethod
    def ThresholdFilter(species, **kwargs):
        """
        Filter used for range-based selection of particles at the current timestep.
        
        Paramters:
        
            species : The particle species object to be filtered (set automatically
                      by the  WarpVisItFilteredSpecies object)
            **kwargs : Define thresholds to be applied via a series of
                       keyword arguments defined as follows var__op=val
                       where var is the variable name, op is the selection
                       operation, and val is the selection value. All 
                       selections are combined via AND, i.e., only elements
                       that suffice all conditions pass the filter. Valid 
                       operations are:
                       
                       * greater
                       * greater_equal,
                       * less
                       * less_equal
                       * equal
                       * not_equal
                       
                       Valid variables include all variables exposed by the
                       WarpVisItFilteredSpecies object via corresponding get
                       function. 
        
        Returns:
        
            indexlist : numpy array of the selected array indicies 
        
        Examples:
        
            * To find all particle with px>1e9 and x<10 define:
              ThresholdFilter(species,
                              px__greater=1e9,
                              x__less=10)
                       
        """
        #Return all indicies if no selection has been applied
        if len(kwargs) == 0:
            numpart = species.getx(gather=0).shape[0]
            return numpy.arange(numpart)
        #Iterate over all selection parameters  and compute the selection
        index = 0
        for key, value in kwargs.items():
            # Get the variable and operation and associated functions
            variable, operation = key.split("__")
            getname = "get"+variable
            variabledata = getattr(species, getname)(gather=0)
            operationfunc = getattr(numpy, operation)
            if index == 0 :
                selection = operationfunc(variabledata, value)
            else:
                selection &= operationfunc(variabledata, value)
            index += 1
            
        #Get the indicies of the selected values
        return numpy.flatnonzero(selection)

    
    #----------------------------------------------------------------------------
    @staticmethod
    def AccumulativeThresholdFilter(self, species, **kwargs):
        """
        This filter is similar to the ThresholdFilter with the main difference that
        the selection is accumulative, i.e., a particle is selected if it suffices
        the threshold condition at the current timestep or did so at any of the 
        previous timesteps at which the filter was evaluated. This function uses 
        the ThresholdSelection filter function and is very similar in it use.
        
        Paramters:
        
            self : With the self paramater as first argument, the filter 
                   function will be bound to the  WarpVisItFilteredSpecies object
                   it is assigned to. This allows the filter function to save 
                   the partilcls IDs selected at previous timesteps. The self
                   paramters is set automatically as usual for bound methods.
                   (Set automatically by the  WarpVisItFilteredSpecies object)
            species : The particle species object to be filtered (set automatically
                      by the  WarpVisItFilteredSpecies object)
            **kwargs : Define thresholds to be applied via a series of
                       keyword arguments defined as follows var__op=val
                       where var is the variable name, op is the selection
                       operation, and val is the selection value. All 
                       selections are combined via AND, i.e., only elements
                       that suffice all conditions pass the filter. Valid 
                       operations are:
                       
                       * greater
                       * greater_equal,
                       * less
                       * less_equal
                       * equal
                       * not_equal
                       
                       Valid variables include all variables exposed by the
                       WarpVisItFilteredSpecies object via corresponding get
                       function. 
        
        Returns:
        
            indexlist : numpy array of the selected array indicies
        
        Examples:
        
            * To find all particle that exceeded a momentum pz>=1e10 at any 
              point in time and that remained within a given transverse
              range of -5<x<5 we can sepcify :
              AccumulativeThresholdFilter(species,
                                          pz__greater_equal=1e10,
                                          x__less_equal=5,
                                          x__greater=-5)
                       
        """
    
        # Compute the list of particls ids selected at the current timestep
        indexlist = WarpVisItSpeciesFilters.ThresholdFilter(species, **kwargs)
        speciespid = species.getpid(gather=0)
        numpart = speciespid.shape[0]
        currentids = speciespid[indexlist]
        # Get the array of previously selected particle ids
        try:
            previousids = numpy.asarray([])
        except AttributeError:
            previousids = self.__FilterIds
        print (currentids, previousids)
        # Merge the current and previous list of particle ids
        self.__FilterIds =  numpy.unique(numpy.hstack((currentids,previousids)),
                                         return_index=False,
                                         return_inverse=False)
        numids = self.__FilterIds.shape[0]
        # Convert the list of ids to a index selection for the current timestep
        return numpy.flatnonzero(numpy.in1d(speciespid, self.__FilterIds))
        
    
    #----------------------------------------------------------------------------
    @staticmethod
    def PidFilter(self, species, pidlist):
        """
        Select particles based on a user-defined list of particle particle ids (pids).
        
        Paramters:
        
            species : The particle species object to be filtered (set automatically
                      by the  WarpVisItFilteredSpecies object)
            pidlist : List of user-defined particle ids (pids) to be selected
        
        Returns:
        
            indexlist : numpy array of the selected array indicies
        """
        return numpy.flatnonzero(numpy.in1d(species.getpid(gather=0), pidlist))
    
    

    #----------------------------------------------------------------------------
    @staticmethod
    def HaloFilter(species, filterRange=3.):
        """
        Filter commonly used to remove halo particles from a particle beam.
        Retrieve particles that are within filterRange*std range of
        the average in x,y,z,ux,uy,uz space.
        
        Paramters:
        
            species : The particle species object to be filtered (set automatically
                      by the  WarpVisItFilteredSpecies object)
            filterRange : Multiple of standard deviation filter range range
            
        Returns:
            
            List of inidices of particles to be selected
        """
        from parallel import globalvar,globalave
        from numpy import sqrt,shape,take,arange,compress,abs
        import sys
        # get particle data
        x = species.getx(gather=0)
        y = species.gety(gather=0)
        z = species.getz(gather=0)
        ux = species.getux(gather=0)
        uy = species.getuy(gather=0)
        uz = species.getuz(gather=0)
        # get rms values
        xrms = sqrt(globalvar(x))
        yrms = sqrt(globalvar(y))
        zrms = sqrt(globalvar(z))
        uxrms = sqrt(globalvar(ux))
        uyrms = sqrt(globalvar(uy))
        uzrms = sqrt(globalvar(uz))
        # get average values
        xave=globalave(x)
        yave=globalave(y)
        zave=globalave(z)
        uxave=globalave(ux)
        uyave=globalave(uy)
        uzave=globalave(uz)
        n=shape(x)[0]
        #sys.stderr.write('rms = %g %g %g %g %g %g\nave = %g %g %g %g %g %g\nn=%d\n'%(xrms,yrms,zrms,uxrms,uyrms,uzrms,xave,yave,zave,uxave,uyave,uzave,n))
        # get indices of particles inside 3*rms in 6-D phase-space
        ii=compress((abs(x-xave)<filterRange*xrms) & \
                    (abs(y-yave)<filterRange*yrms) & \
                    (abs(z-zave)<filterRange*zrms) & \
                    (abs(ux-uxave)<filterRange*uxrms) & \
                    (abs(uy-uyave)<filterRange*uyrms) & \
                    (abs(uz-uzave)<filterRange*uzrms) \
                    ,arange(n))
        #sys.stderr.write('%s\n'%(str(ii)))
        return ii






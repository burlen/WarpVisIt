import sys
import os
import parallel
import time

#-----------------------------------------------------------------------------
def pStatus(msg):
    """
    Rank 0 prints the message to sderr
    """
    if parallel.get_rank() == 0:
        tstr = time.strftime("%H:%M:%S", time.localtime())
        sys.stderr.write('=====: %d %s %s \n'%(parallel.get_rank(), tstr, msg))
        sys.stderr.flush()
    return

#-----------------------------------------------------------------------------
def pDebug(msg):
    """
    Print an debug message to stderr (if the following line is not commented)
    """
    #tstr = time.strftime("%H:%M:%S", time.localtime())
    #sys.stderr.write('=====: %d %s %s \n'%(parallel.get_rank(), tstr, msg))
    #sys.stderr.flush()
    return

#-----------------------------------------------------------------------------
def pError(msg):
    """print an error message to stderr"""
    tstr = time.strftime("%H:%M:%S", time.localtime())
    sys.stderr.write('Error: %d %s %s\n'%(parallel.get_rank(), tstr, msg))
    sys.stderr.flush()
    return

#-----------------------------------------------------------------------------
def getEnvVar(name, retType=str, defaultVal=None):
    """
    Get the environmant variable and cast to the desired
    type otherwise return the default value
    """
    var = os.getenv(name)
    if var:
        return retType(var)
    else:
        return defaultVal

##############################################################################
class VisItEnv:
    """
    Scrape the runtime environment for VISITARCH and
    build the neccessary paths.
    """
    __sharedState = {}

    #-------------------------------------------------------------------------
    def __init__(self):
        self.__dict = self.__sharedState
        if not self.__sharedState:
            self.__root = os.getenv('VISIT')
            if not self.__root:
                raise RuntimeError('VISIT is not set')
            self.__bin = os.path.join(self.__root, 'bin')
            self.__lib = os.path.join(self.__root, 'lib')
            self.__sitePackages = os.path.join(self.__lib, 'site-packages')
            self.__Validate()
            sys.path.append(self.__sitePackages)
            sys.path.append(self.__lib)

    #-------------------------------------------------------------------------
    def __Dict(self):
        """Make a dictionary from the internal state"""
        return {'root':self.__root,
            'bin':self.__bin,
            'lib':self.__lib,
            'site-packages':self.__sitePackages
            }

    #-------------------------------------------------------------------------
    def __repr__(self):
        """Convert to human readable string"""
        s = ''
        for name,path in self.__Dict().iteritems():
            s += '%s = %s\n'%(name,path)
        return s

    #-------------------------------------------------------------------------
    def __Validate(self):
        """Check paths and raise an error if any are bad"""
        for name,path in self.__Dict().iteritems():
            if not os.path.exists(path):
                raise RuntimeError('Enviroment is not configured. %s=%s  doesn\'t exist'%(name,path))

    #-------------------------------------------------------------------------
    def GetRoot(self):
        """Return the VisIt's __root dir"""
        return self.__root

    #-------------------------------------------------------------------------
    def GetLib(self):
        """Return the VisIt's library dir"""
        return self.__lib

    #-------------------------------------------------------------------------
    def GetBin(self):
        """Return the VisIt's binrary dir"""
        return self.__bin

    #-------------------------------------------------------------------------
    def GetSitePackages(self):
        """Return the VisIt's python site-packages dir"""
        return self.__sitePackages

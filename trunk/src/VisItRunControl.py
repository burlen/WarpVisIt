
class VisItRunControl:
    """
    Some data used to control run
    """
    MODE_START = 0
    MODE_RUN = 1
    MODE_STOP = 2
    def __init__(self):
        """ """
        self.__NumIteration = 1
        self.__Iteration = 0
        self.__VisInterval = 2
        self.__SimTime = 0
        self.__Mode = MODE_START

    def Step(self):
        """ """
        self.__Iteration += 1

    def SetSimulationTime(self, time):
        """ """
        self.__SimTime = time

    def SetVisualizationInterval(self, num):
        """ """
        self.__VisInterval = num

    def SetNumberOfIterations(self, num):
        """ """
        self.__NumIteration = num

    def Visualize(self):
        """ """
        return self.__Iteration % self.__VisInterval == 0

    def RunComplete(self):
        """ """
        return self.__Iteration >= self.__NumIteration

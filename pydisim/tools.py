
class Recorder:
    '''
    This class represents a collection of temporal data that is recorded at
    specific timestamps.
    '''
    def __init__(self,names):
        '''
        The data structure is initialised with a list of variable names. The
        names will be used to decorate plots and can be anything. The size of
        the data structure is determined from the length of names

        Arguments:
        ----------
        names : list of strings
            List of nice names of variables.
        '''

        self.names = names
        self.data = []
        for n in names:
            self.data.append([])


    def record(self,*args):
        '''
        Record a snapshot

        Arguments:
        ----------
        *args : list
            All the values to record. the number of values must be equal to the
            size of the structure. Additional data will be ignored.
        '''

        for i in range(len(self.data)):
            if len(args) <= i:
                self.data[i].append(None)
            else:
                self.data[i].append(args[i])

    def plot(self,ax,timestamps=None):
        '''
        Plot recorder's data on matplotlib axis ax

        Arguments:
        ----------
        ax : matplotlib axis
            Matlpotlib axis to plot data on
        timestamps : list
            List of x-axis values
        '''

        if timestamps == None:
            for i in range(len(self.data)):
                ax.plot(self.data[i],label=self.names[i])
        else:
            for i in range(len(self.data)):
                ax.plot(timestamps,self.data[i],label=self.names[i])
        ax.legend()
        

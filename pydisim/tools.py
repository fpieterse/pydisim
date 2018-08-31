import matplotlib.pyplot as plt

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

class PIDRecorder():
    '''
    A special type of recorder used to record a PID controller. It doesn't need
    variable names when instantiated (they are always SP,PV,OP). It also makes
    it easier to plot the SP/PV on the same axis and the OP on a unique axis
    '''

    def __init__(self,tagname=None):
        '''
        Arguments:
        ----------
        tagname : string
            The tagname of the PID controller
        '''
        self.sppvRec = Recorder(['SP','PV'])
        self.opRec = Recorder(['OP'])
        self.tagname = tagname

    def record(self,sp,pv,op):
        self.sppvRec.record(sp,pv)
        self.opRec.record(op)

    def plot(self,ax1,ax2,timestamps=None):
        '''
        Plot recorders data on matplotlib axes

        Arguments:
        ----------
        ax1, ax2 : matplotlib axis
            SP and PV will be plotted on ax1 and OP will be plotted on ax2. If
            these are set to None then a new figure will be created.
        timestamps : list
            List of x-axis values
        '''
        if ax1 == None:
            fig,axes = plt.subplots(2,1,sharex=True)
            ax1 = axes[0]
            ax2 = axes[1]
            fig.tight_layout()
        else:
            fig = None

        self.sppvRec.plot(ax1,timestamps)
        self.opRec.plot(ax2,timestamps)
        ax1.set_ylabel(self.tagname)


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

        if type(names) != list:
            raise ValueError('Expected list of variable names')

        self.names = names
        self.clear()


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

    def plot(self,ax=None,timestamps=None):
        '''
        Plot recorder's data on matplotlib axis ax

        Arguments:
        ----------
        ax : matplotlib axis
            Matlpotlib axis to plot data on
        timestamps : list
            List of x-axis values
        '''

        if ax == None:
            fig = plt.figure()
            ax = plt
        if timestamps == None:
            for i in range(len(self.data)):
                ax.plot(self.data[i],label=self.names[i])
        else:
            for i in range(len(self.data)):
                ax.plot(timestamps,self.data[i],label=self.names[i])
        ax.legend()

    def clear(self,keep=None):
        '''
        Clear the recorded history

        Arguments:
        ----------
        keep : int
            Number of samples to keep. If None the whole thing is cleared
        '''
        if keep is not None and (keep > 0):
            for (i,d) in enumerate(self.data):
                # this does not raise OutOfRangeExceptionata
                self.data[i] = d[-keep:] 
        else:
            self.data = []
            for i in range(len(self.names)):
                self.data.append([])



class PIDRecorder(Recorder):
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
        super().__init__([tagname + '.sp',
                          tagname + '.pv',
                          tagname + '.op'])
        self.tagname = tagname

    def record(self,pid):
        super().record(pid.sp,pid.pv,pid.op)

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

        if timestamps is None:
            for (i,ax) in zip((0,1,2),(ax1,ax1,ax2)):
                ax.plot(self.data[i],label=self.names[i])
        else:
            for (i,ax) in zip((0,1,2),(ax1,ax1,ax2)):
                ax.plot(timestamps,self.data[i],label=self.names[i])

        ax1.set_ylabel(self.tagname)


def save_to_file(filename,time,*recorders):
    '''
    Save recorder data to file.

    Arguments:
    ----------
    filename : string
        name of file
    time : list
        list of timestamps corresponding to recorder data
    *recorders : Recorder
        Remaining arguments are recorders to be saved to file
    '''

    with open(filename,'w') as f:
        f.write('time')
        for r in recorders:
            for n in r.names:
                f.write(',' + n)
        f.write('\n')

        for i in range(len(time)):
            f.write(str(time[i]))
            for r in recorders:
                for d in r.data:
                    f.write(',{:.6f}'.format(d[i]))
            f.write('\n')

                

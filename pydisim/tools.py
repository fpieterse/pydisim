import matplotlib.pyplot as plt
import numpy
import datetime
import pandas

from .processes import AbstractProcess


class RecorderProcess(AbstractProcess):
    # All recorders have the same starting time, set when the first one is
    # created
    def __init__(self,names,rec_int=0):
        '''
        Arguments:
        ----------
        names : list of strings
            Names of inputs.  A parameter will be created on the class with this
            name so you have to make sure it does not overwrite any existing
            values.
        rec_int : int
            Minimum recording interval (seconds).
        max_hist : float
            Maximum recording duration (hours)
        '''
        super().__init__()
        self.process_manager._recorders.append(self)
        

        self.names = names
        self.rec_int = rec_int

        for name in names:
            try:
                attr = getattr(self,name)
                # if no error, we found the attribute
                # TODO: Write test for this
                raise Exception("Attribute " + name + " already exists")
            except AttributeError:
                setattr(self,name,numpy.nan)
        self.last_t = 0

        self.data = {name:[] for name in names}
        self.timestamps = []


    def run_for(self,dt):
        # Recorder ignores dt, check process_manager to see where it is at
        # This is to keep different recorder syncronised.

        t = self.process_manager.simulation_time
        if self.last_t + self.rec_int <= t:
            self.timestamps.append(self.process_manager.simulation_start +
                                   datetime.timedelta(seconds = t))
            for name in self.data:
                self.data[name].append( getattr(self,name) )

            self.last_t = t

    def plot(self,ax=None,subplots=False):
        '''
        Plot recorder's data on matplotlib axis ax

        Arguments:
        ----------
        ax : matplotlib axis
            Matlpotlib axis to plot data on
        subplots : bool
            Set to True to plot values as subplots. Only used if ax==None

        '''

        fig = None
        if ax == None:
            if subplots == True:
                nr = len(self.data)
                fig,ax = plt.subplots(nrows=nr,sharex=True)
            else:
                fig,ax = plt.subplots(1,1)
        else:
            subplots = False
       
        i = 0
        for key in self.data:
            if subplots:
                _ax = ax[i]
            else:
                _ax = ax
            _ax.plot(self.timestamps,self.data[key],label=key)
            if subplots:
                _ax.set_ylabel(key)
            i += 1

        if not subplots:
            ax.legend()

        if fig != None:
            fig.tight_layout

        

    def clear(self,keep=None):
        '''
        Clear the recorded history

        Arguments:
        ----------
        keep : int
            Number of minutes to keep. If None the whole thing is cleared
        '''
        if keep is not None and (keep > 0):
            if self.rec_int <= 0:
                self.rec_int = 1
            n_keep = int(keep * self.rec_int / 60)
            for key in self.data:
                # this does not raise OutOfRangeExceptionata
                del self.data[key][:-n_keep]
                del self.timestamps[:-n_keep]
        else:
            for key in self.data:
                self.data[key].clear()
                self.timestamps.clear()

    def as_dataframe(self):
        '''
        Return data as pandas dataframe
        '''
        df = pandas.DataFrame(
            index = self.timestamps,
            data = self.data)
        return df

        


class PIDRecorderProcess(RecorderProcess):
    def __init__(self,loopname,PIDprocess,dt):
        self.loopname=loopname
        self.pid = PIDprocess
        names = [ loopname+suf for suf in ['.sp','.pv','.op'] ]
        super().__init__(names,dt)


    def run_for(self,dt):
        for par in ['sp','pv','op']:
            name = self.loopname + '.' + par
            setattr(self,name,getattr(self.pid,par))

        super().run_for(dt)

    def plot(self,ax=None):
        if ax == None:
            fig,ax = plt.subplots(1,1)

        ax2 = ax.twinx()

        lines = []
       
        lines += ax.plot(self.timestamps,
                     self.data[self.loopname + '.pv'],
                     label=self.loopname + '.pv')
        lines += ax.plot(self.timestamps,
                     self.data[self.loopname + '.sp'],
                     label=self.loopname + '.sp')
        lines += ax2.plot(self.timestamps,
                     self.data[self.loopname + '.op'],
                     label=self.loopname + '.op',
                     c='C2')

        labels = [line.get_label() for line in lines]
        ax.legend(lines,labels)


        ax.set_title(self.loopname)



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

    recList = []
    for r in recorders:
        if isinstance(r,Recorder):
            recList.append(r)
        elif type(r) == dict:
            for key in r:
                if isinstance(r[key],Recorder):
                    recList.append(r[key])
                else:
                    print('WARNING: Element in dictionary is not a Recorder')
                    print('         Data will not be saved to file')
        elif type(r) == list:
            for ri in r:
                if isinstance(ri,Recorder):
                    recList.append(ri)
                else:
                    print('WARNING: Element in list is not a Recorder')
                    print('         Data will not be saved to file')
        else:
            print('WARNING: Argument is not a Recorder')
            print('         Data will not be saved to file')

    with open(filename,'w') as f:
        f.write('time')
        for r in recList:
            for n in r.names:
                f.write(',' + n)
        f.write('\n')

        for i in range(len(time)):
            f.write(str(time[i]))
            for r in recList:
                for d in r.data:
                    f.write(',{:.6f}'.format(d[i]))
            f.write('\n')

                

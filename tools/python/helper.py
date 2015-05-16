import numpy as np
import pickle

goldenratio=0.61803398875

def dash(index):
    ''' centrally controlls dashes for plots'''
    mydashes = [(None,None),(1,3),(2,4,7,4),(6,3),(2,4,2,4,7,4),(7,4,7,4,2,4),
                (2,4,2,4,2,4,7,4),(2,4,2,4,7,4,7,4)]
    nelem = len(mydashes)
    
    return mydashes[np.mod(index,nelem)]


def tickLabelFontSize(ax,fontsize):

        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
        ax.yaxis.get_offset_text().set_fontsize(fontsize-2)


def tickLabelFontSizeCb(cb,fontsize,rasterized=True):
    for tick in cb.ax.get_yticklabels():
        tick.set_fontsize(fontsize)
    cb.ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    cb.ax.yaxis.get_offset_text().set_fontsize(fontsize-2)
    cb.solids.set_rasterized(rasterized) 

def scale_to_unit_interval(ndar):
    """ Scales all values in the ndarray ndar to be between 0 and 1 """
    ndar = ndar.copy()
    ndar -= ndar.min()
    ndar *= 1.0 / (ndar.max())
    return ndar


def dump(data,dumpfile='dumpfile.dat'):

    file=open(dumpfile,'wb') 
    pickle.dump(data,file, pickle.HIGHEST_PROTOCOL)
    file.close()


def restore(dumpfile='dumpfile.dat'):

    file=open(dumpfile,'rb')
    data = pickle.load(file)
    file.close()
    return data

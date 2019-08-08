#!/usr/bin/python
# -*- coding: UTF-8 -*-
from __future__ import absolute_import
from __future__ import with_statement
from __future__ import division
from __future__ import nested_scopes
from __future__ import generators
from __future__ import unicode_literals
from __future__ import print_function

import os,sys,traceback,h5py,gc

# Define paths to data files and disk cache
path      = os.path.abspath(os.path.expanduser('~/Workspace2/PPC_data/'))+os.sep
CACHEDIR  = os.path.abspath(os.path.expanduser('./'))
CACHENAME = 'cache'
print('Data  location is',path)
print('Cache location is %s/%s'%(CACHEDIR,CACHENAME))

# Disk cache some intermediate results 
from neurotools.jobs.initialize_system_cache import initialize_caches
initialize_caches(
    level1=CACHEDIR,
    force  =False,
    verbose=False,
    CACHE_IDENTIFIER = CACHENAME)
import neurotools
memoize = neurotools.jobs.ndecorator.memoize

# helper functions 
from neurotools.nlab   import *
from neurotools.hdfmat import printmatHDF5, hdf2dict, getHDF

# Linear system solvers
from scipy.linalg import solve,lstsq

# Get Numpy and matplotlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Get helper functions
# from ppc_util import *
# import ppc_trial

# override depricated matplotlib `find` and `amap` functions
from  neurotools.tools import find, amap

def get_file(animal,session):
    '''
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    filename : str
        Data file name as a string
    '''
    return 'm%02d_s%02d.mat'%(animal,session)

def get_data(animal,session):
    '''
    Get HDF5 file instance corresponding to given animal and session.

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return h5py.File(path+get_file(animal,session), 'r')

def release_files(clear_cache=False):
    '''
    Ad-hoc way to deal with poor memory management. 
    
    This clears some of the memoization caches and also forces
    garbage collections and closes any files that might have 
    accidentally been left open.
    
    Other Parameters
    ----------------
    clear_cache : bool, default `False`
        Whether to clear the RAM caches as well as the open files.
    '''
    # Browse through ALL objects
    for obj in gc.get_objects():   
        try:
            # Close any open HDF5 files
            if isinstance(obj, h5py.File):   
                try:    obj.close()
                except: pass # Was already closed
        except: pass
    if clear_cache: clear_memoized()
    gc.collect()

def get_dFF(animal,session,units=None):
    '''
    Extract the 'dFF' variable from a PPC experiment
    
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Other Parameters
    ----------------
    units : integer array-like
        List of units for which to return the dFF traces. Defaults to 
        returning all units (including ones with missing data). 
    
    Returns
    -------
    dFF : Ntimes x Nneurons array
        Raw dF/F calcium signals. First axis is time-points and second
        axis is neurons.
    '''
    hdfdata = get_data(animal,session)
    neural  = getHDF(hdfdata,'session_obj/timeSeries/calcium/data')
    dFF,decon_Ca = neural
    hdfdata.close()
    return dFF if units is None else dFF[:,np.int32(units)]

@memoize
def get_smoothed_dFF(animal,session,unit,Fmin,Fmax,zeromean=True):
    '''
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
    unit : int
        Which (0-indexed) unit ID to use
    Fmin : float
        Low-frequency cutoff for filtering
    Fmax : float
        High-frequency cutoff for filtering

    Returns
    -------
    dFF : Ntimes x Nneurons array
        Filtered dF/F calcium signals. 
        First axis is time-points and second axis is neurons.
    '''
    z  = get_dFF(animal,session,unit)
    FS = get_FS(animal,session)
    z  = bandpass_filter(z,fa=Fmin,fb=Fmax,Fs=FS)
    if zeromean:
        z = z-np.mean(z)
    return z

@memoize
def get_subject_ids():
    '''
    Report a list of subjects that have data

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    #global path
    fs = os.listdir(path)
    alli = []
    for m in os.listdir(path):
        m = m.lower()
        if not m.endswith('mat'): continue
        m = m.split('.mat')[0]
        if not '_' in m: continue
        m = m.split('_')[0]
        try:
            alli.append(int(m[1:]))
        except:
            pass
    alli = sorted(list(np.unique(alli)))
    return alli

@memoize
def get_session_ids(animal):
    '''
    Report a list of available sesson numbers for a given subject.
    
    Some invalid sessions have been hard-coded here to be removed.

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    sessions : list
        List of available session IDs
    '''
    prefix = 'm%02d_'%animal
    fs = os.listdir(path)
    ss = []
    for m in [f for f in fs if f.startswith(prefix)]:
        m = m.split('.mat')[0].split(prefix)[1]
        if m.startswith('s') and len(m)==3:
            s = int(m[1:])
            if 'm%02d_s%02d.mat'%(animal,s) in os.listdir(path):
                ss.append(s)
    return sorted(list(set(ss)))

@memoize
def good_units(animal,session):
    '''
    Rather than use the confidence label stored with the dataset,
    this routine actually checks which units have available data.
    This is pretty slow, so results are cached on disk.

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    good : np.array
        Boolean list; 
        available: `True`, unavailable: `False`
    '''
    y    = get_dFF(animal,session)
    good = np.all(np.isfinite(y),axis=0)
    return good

def good_units_index(animal,session):
    '''
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    list of indecies for units available on the given session
    '''
    return find(good_units(animal,session))

def get_FS(animal,session):
    '''
    Get sampling rate of given animal and session.

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    FS : float
        The sampling frame rate of the session
    '''
    hdfdata = get_data(animal,session)
    FS = getHDF(hdfdata,'session_obj/timeSeries/frameRate') # Hz
    hdfdata.close()
    return FS

def get_nonneural(animal,session):
    '''
    Get the non-neural timeseries (kinematics, mostly)
    
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    hdfdata = get_data(animal,session)
    result = getHDF(hdfdata,'session_obj/timeSeries/virmen/data')
    hdfdata.close()
    return result

def get_intrial(animal,session):
    '''
    Gets whether each time-point is in a trial or not. 
    This is based on the "in-trial" flag provided in the recording data.
    
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    np.array
    '''
    return ~(get_nonneural(animal,session).T[8]>0)

def get_timestamp(animal,session):
    '''
    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return get_nonneural(animal,session).T[0]

def get_x(animal,session):
    '''
    x location in virtual maze
    units are meters

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return get_nonneural(animal,session).T[1]

def get_y(animal,session):
    '''
    y location in virtual maze
    units are meters

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return get_nonneural(animal,session).T[2]

def get_theta(animal,session):
    '''
    Get head-direction view angle. Units are in degrees. 
    
    The Virmen system will wrap the angle if the animal is doing 
    corkscrews. This routine add circular wrapping.

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return (get_nonneural(animal,session).T[3]+180)%360-180

def get_theta_radians(animal,session):
    '''
    Head direction in radians

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return (get_theta(animal,session)*np.pi/180+3*pi)%(2*pi)-pi

def get_dx(animal,session):
    '''

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    '''
    return get_nonneural(animal,session).T[4]

def get_dy(animal,session):
    '''
    y velocity. Units are meters per second

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    np.array
    '''
    return get_nonneural(animal,session).T[5]

kininfo = {
    0:{'name':'X position'    ,'units':'m'  ,'get':get_x},
    1:{'name':'Y position'    ,'units':'m'  ,'get':get_y},
    2:{'name':'X velocity'    ,'units':'m/s','get':get_dx},
    3:{'name':'Y velocity'    ,'units':'m/s','get':get_dy},
    4:{'name':'Head direction','units':'°'  ,'get':get_theta}
}

def get_type(animal,session):
    '''
    Get trial type (black right = 2, white left = 3) for given
    animal and session.

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    type : np.array
        Array of trial-type indicators.
        (black right = 2, white left = 3)
    '''
    return get_nonneural(animal,session).T[6]

def get_reward(animal,session):
    '''
    Timeseries indicating reward
    0 = no reward
    1 = reward timepoint

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : int
        Which session ID to use
        
    Returns
    -------
    np.array
    '''
    return get_nonneural(animal,session).T[7]

@memoize
def get_units_in_common(animal,sessions):
    '''
    Get neurons in common to a set of sessions

    Parameters
    ----------
    animal : int
        Which subject ID to use
    session : list of ints
        Which session IDs to use


    Returns
    -------
    units: 
        List of unit numbers in common
    uidxs:
        Index into list of good units for each session
        
    '''
    unitmap  = {s:set(good_units_index(animal,s))      for s in sessions}
    units    = sorted(list(set.intersection(*unitmap.values())))
    uidxs    = [find([u in units for u in unitmap[s]]) for s in sessions]
    return units,uidxs


def get_consecutive_recordings(animal,
                               MINNEURONS=200,
                               MINDURDAYS=4,
                               verbose=True):
    '''
    Identify consecutive spans of days for the given `animal`,
    with at least `MINNEURONS` in common and lasting at least 
    `MINDURDAYS` long.
    
    Parameters
    ----------
    animal: int
        which subject to use
        
    Other Parameters
    ----------------
    MINNEURONS: int
        Minimum number of neurons in common.
    MINDURDAYS: int
        Minimum duration in days.
    verbose: bool
        Print debugging info.
        
    Returns
    -------
    valid_spans : 
        dictionary mapping tuples of starting and ending session 
        numbers to the list of neurons that they share in common.
        
    '''
    # Get session information
    sessions = np.array(get_session_ids(animal))
    if verbose:
        print('Testing subject %d'%animal)
        print('  Available sessions',' '.join(map(str,sessions)))
    # Find all spans of sessions over consecutive days
    d = np.diff(sessions)
    L = len(sessions)
    spans = [(i,j) for i in range(L-1)
             for j in range(i+1,L)
             if all(d[i:j]==1)]
    spans = sessions[np.int32(spans)]
    # Get good units for all included sessions
    unitmap = {s:set(good_units_index(animal,s)) for s in sessions}
    # Get overlapping unit sets for all spans of sessions
    common = {(a,b):set.intersection(*[unitmap[s]
                                       for s in range(a,b+1)])
              for (a,b) in spans}
    # Get durations of valid spans of days
    durations = np.diff(spans).ravel()+1
    if verbose:
        print('  There are %2d spans at least %d days long; of these…'\
              %(sum(durations>=MINDURDAYS),MINDURDAYS))
    spans = spans[durations>=MINDURDAYS,:]
    # Get spans of days with enough units
    valid_spans = np.array([span for span in spans
                         if len(common[tuple(span)])>=MINNEURONS])
    # Hack to remove subsets
    daysets = [set(range(a,b+1)) for (a,b) in valid_spans]
    useme = []
    N = len(valid_spans)
    for i in range(N):
        if not any([daysets[i].issubset(daysets[j]) 
                    for j in range(N) if i!=j]):
            useme.append(valid_spans[i])
    valid_spans = useme
    if verbose:
        print('  There are %d spans with at least %d neurons:'\
              %(len(valid_spans),MINNEURONS))
        [print('  \t%02d-%02d (%d days, %3d neurons)'\
               %(a,b,b-a+1,len(common[a,b]))) 
         for a,b in valid_spans]
    return {k:sorted(list(common[k])) for k in map(tuple,valid_spans)}


def add_constant(data):
    return cat([data,ones((data.shape[0],1))],axis=1)

def polar_error_degrees(x,xh):
    # Report error in physical units
    e = abs(x-xh)
    e[e>180] = 360-e[e>180]
    return mean(abs(e)**2)**.5, mean(abs(e))


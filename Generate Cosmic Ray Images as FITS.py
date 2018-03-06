
# coding: utf-8

# # Let's create a batch of images with cosmic rays
# This notebook makes a bundle of 2-second images that have cosmic rays injected into them. To make it easier to play with, it's made as a subarray.

# In[1]:


# this sets where to look for inputs and store outputs 
# (it will default to `~/.tess/spyffi` if the $SPYFFIDATA environment variable isn't set)
import os
os.environ["SPYFFIDATA"] = '/Volumes/dandelion/Cosmos/Data/TESS/FFIs'


# In[2]:


# this makes sure that updates are actually printed to the terminal (instead of just to a log)
import logging
logging.basicConfig(level="INFO")


# In[3]:


# this loads a default observation definition, which is a dictionary of dictionaries
from SPyFFI.Observation import Observation, default
inputs = default


# In[5]:


# give a label to this observation (that goes into the directory name)
inputs['camera']['label'] = 'manycosmics'
# how big of a subarray in pixels (centered on the FOV); None will give four normal CCDs
inputs['camera']['subarray'] = 400
# should we change the focus throughout the orbit? 
inputs['camera']['variablefocus'] = False

# should we skip injecting cosmic rays? (False = *do* inject cosmics)
inputs['expose']['skipcosmics'] = False
# let's write the cosmic ray images out as separate files
inputs['expose']['writecosmics'] = True
inputs['expose']['writenoiseless'] = False
inputs['expose']['writesimulated'] = True
inputs['expose']['compress'] = {2: False, 20: False, 120: False, 1800: False}

# what's the catalog ("testpattern" or "UCAC4")
inputs['catalog']['name'] = 'testpattern'
# should we randomize the magnitudes of the stars
inputs['catalog']['testpatternkw']['randomizemagnitudes'] = True
# what range of magnitudes should the stars span?
inputs['catalog']['testpatternkw']['magnitudes'] = [6,16]

# how many of each exposure duration should we make?
days = 1.0
N = int(days*24*60*60/2)
inputs['observation']['cadencestodo'] = {2:N}
print("Going to make {}".format(inputs['observation']['cadencestodo']))


# In[6]:


o = Observation(inputs)


# In[7]:


# generate some light curves to populate the catalog
import SPyFFI.Lightcurve as lc, numpy as np

# every star in the catalog has a light curve; we can modify them one-by-one
for i in range(len(o.camera.catalog.ra)):
    o.camera.catalog.lightcurves[i] = lc.constant()

# we can convert the catalog RA + Dec to x, y
ra, dec = o.camera.catalog.ra, o.camera.catalog.dec
coord = o.camera.cartographer.point(ra, dec, 'celestial')
x, y = coord.ccdxy.tuple

# let's pick one bright star that's on the detector...
n = o.camera.ccds[0].xsize
ok = (x > 0)*(x < n)*(y < n)*(y > 0)
print(sum(ok))


# In[8]:


# let's print out the locations and light curves of the stars
for i, l in enumerate(o.camera.catalog.lightcurves):
    print '{:>5} = ({:5.1f},{:5.1f}) = {}'.format(i,x[i],y[i],l)


# In[9]:


o.create()


# This should now take a while (like a few minutes?) to create a bunch of images and drop them in `$SPYFFIDATA/testpattern_6to16_transit/1800s/sub200x200/`. Go check them out and see what's what!

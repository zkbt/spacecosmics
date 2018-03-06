
# coding: utf-8

# # Creating Subarray Images with SPyFFI
# This notebook shows some of the basics of SPyFFI, using it to create images over a small subarray of the detector. Because this is rendering fewer pixels, it's probably the fastest way to test `SPyFFI` end-to-end to make sure it's working.

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


# The inputs are all organized into dictionaryies within the `inputs` dictionary.

# In[4]:


# how do we define the catalog of stars to use (and the light curves)
print default['catalog']


# In[5]:


# how should exposures be created (what should be written out?)
print default['expose']


# In[6]:


# set up the geometry of the camera, and (importantly) the PSF library
print default['camera']


# In[7]:


# in what ways should we jitter the motion of the camera?
print default['jitter']


# In[8]:


# how many of each exposure time should we create (and in what order?)
print default['observation']


# Let's modify some of those input options, by changing values associated with some of the directory keys inside each sub-group.

# In[9]:


# give a label to this observation (that goes into the directory name)
inputs['camera']['label'] = 'manycosmics'
# how big of a subarray in pixels (centered on the FOV); None will give four normal CCDs
inputs['camera']['subarray'] = 200
# should we change the focus throughout the orbit? 
inputs['camera']['variablefocus'] = False

# should we skip injecting cosmic rays? (False = *do* inject cosmics)
inputs['expose']['skipcosmics'] = False
# let's write the cosmic ray images out as separate files
inputs['expose']['writecosmics'] = True

# what's the catalog ("testpattern" or "UCAC4")
inputs['catalog']['name'] = 'testpattern'
# should we randomize the magnitudes of the stars
inputs['catalog']['testpatternkw']['randomizemagnitudes'] = True
# what range of magnitudes should the stars span?
inputs['catalog']['testpatternkw']['magnitudes'] = [6,16]

# how many of each exposure duration should we make?
inputs['observation']['cadencestodo'] = {2:100}
o = Observation(inputs)


# In[10]:


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
ok *= (o.camera.catalog.tmag < 11.0)
lucky = np.random.choice(np.nonzero(ok)[0], 1)[0]

# ... and inject a deep transit into it
o.camera.catalog.lightcurves[lucky] = lc.Trapezoid(P=0.3, 
                                                   E=2457827.0 + 0.15, 
                                                   D=1.0, 
                                                   T23=0.00, 
                                                   T14=0.1)


# In[11]:


# let's print out the locations and light curves of the stars
for i, l in enumerate(o.camera.catalog.lightcurves):
    print '{:>5} = ({:5.1f},{:5.1f}) = {}'.format(i,x[i],y[i],l)


# In[12]:


# let's plot the light curve of that one lucky star
get_ipython().run_line_magic('matplotlib', 'inline')
l = o.camera.catalog.lightcurves[lucky]
l.demo()
print l


# In[13]:


o.create()


# This should now take a while (like a few minutes?) to create a bunch of images and drop them in `$SPYFFIDATA/testpattern_6to16_transit/1800s/sub200x200/`. Go check them out and see what's what!

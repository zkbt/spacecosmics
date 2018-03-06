'''Generate TESS pixel lightcurve cubes with dimensions (xpix)x(ypix)x(time).'''
from imports import *
#from Strategies import *
from Stacker import Central, Sum
from SPyFFI.Observation import Observation, default
from craftroom.displays.ds9 import ds9


# create a default input dictionary
cosmicinputs = default
# give a label to this observation (that goes into the directory name)
cosmicinputs['camera']['label'] = 'cube'

# should we change the focus throughout the orbit?
cosmicinputs['camera']['variablefocus'] = False
cosmicinputs['jitter']['amplifyinterexposurejitter'] = 1.0

# should we skip injecting cosmic rays? (False = *do* inject cosmics)
cosmicinputs['expose']['skipcosmics'] = False

# let's write the cosmic ray images out as separate files
cosmicinputs['expose']['writecosmics'] = False
cosmicinputs['expose']['writenoiseless'] = False
cosmicinputs['expose']['writesimulated'] = False
cosmicinputs['expose']['compress'] = {2: False, 20: False, 120: False, 1800: False}

# what's the catalog ("testpattern" or "UCAC4")
cosmicinputs['catalog']['name'] = 'testpattern'
# should we randomize the magnitudes of the stars
cosmicinputs['catalog']['testpatternkw']['randomizemagnitudes'] = True
# what range of magnitudes should the stars span?
cosmicinputs['catalog']['testpatternkw']['magnitudes'] = [10,10]




class Cube(Talker):
	'''Cube to handle simulated postage stamp pixel light curves;
			has dimensions of (xpixels, ypixels, time).'''

	def __init__(self, subarray=42, n=15, cadence=120, subexposurecadence = 2, inputs=cosmicinputs, stacker=Central(10), **kwargs):
		'''Initialize a cube object.

		keyword arguments:
		[inputs] = dictionary of inputs to a SPyFFI observation
		[stacker] = a stacker object for combining images together
		[**kwargs] = what happens to these?
		'''

		# decide whether or not this Cube is chatty
		Talker.__init__(self)

		# keep track of the inputs
		self.inputs = inputs

		self.inputs['observation']['cadencestodo'] = {cadence:n}
		self.inputs['camera']['subarray'] = subarray
		assert(subarray is not None)


		for k in ['writesimulated', 'writenoiseless', 'writesimulated']:
			self.inputs['expose'][k] = False

		# keep track of which stacker is being used to create this image
		self.stacker = stacker

		# create an observation
		self.cadence = cadence
		self.observation = Observation(inputs)
		self.camera = self.observation.camera
		self.ccd = self.observation.camera.ccds[0]

		# define basic geometry
		self.size = subarray
		self.xpixels = self.size
		self.ypixels = self.size

		# number of images and cadence
		self.n = n
		self.subexposurecadence = subexposurecadence
		self.shape = (self.xpixels, self.ypixels, self.n)


		# do we need to create a subcube of subexposures?
		self.ninstack = np.round(self.cadence/self.subexposurecadence).astype(np.int)
		if self.ninstack > 1:
			self.subcube = Cube(subarray=self.size,
								n=self.ninstack,
								cadence=subexposurecadence,
								inputs=self.inputs,
								stacker=Sum())
			self.subcube.camera.catalog = self.camera.catalog
			self.camera.setCadence(subexposurecadence)


		# create empty (xpixels, ypixels, n)
		if self.cadence == 2:
			bits = np.float32
		else:
			bits = np.float64

		#
		self.photons, self.cosmics, self.noiseless, self.unmitigated = np.zeros(self.shape).astype(bits), np.zeros(self.shape).astype(bits), np.zeros(self.shape).astype(bits), np.zeros(self.shape).astype(bits)

		# create a dictionary to store a bunch of summaries
		self.summaries = {}

		# populate the cube with simulated pixel data
		# self.simulate()
		# self.plot()

	def __repr__(self):
		return '<Cube | {} | {}s>'.format(self.shape, self.cadence)

	def simulate(self):

		# if we're using a "Sum" stacking strategy, then don't generate the individual 2-second frames
		if (self.stacker.name == 'Sum') | (self.cadence == self.subexposurecadence):
			# loop over (already stacked) exposures, creating them
			for i in range(self.n):
				self.speak('filling exposure #{0:.0f}/{1:.0f}'.format(i, self.n))
				p, c, n = self.ccd.expose(**self.inputs['expose'])
				self.photons[:,:,i], self.cosmics[:,:,i], self.noiseless[:,:,i] = p, c, n
			self.unmitigated = self.photons

			# store some useful accessory information about
			self.background = self.ccd.backgroundimage
			self.noise = self.ccd.noiseimage
			self.catalog = self.camera.catalog
		else:
			# otherwise, we need to create the subcube
			theWayToStack = self.stacker
			self.speak('using {0} strategy to stack from 2-second images')
			for i in range(self.n):
				self.speak()
				self.speak('creating a temporary stack of {0} images {1}s images'.format(self.ninstack, self.subexposurecadence))
				self.subcube.simulate()
				p, c, n, u = theWayToStack.stack(self.subcube, self.ninstack)
				self.photons[:,:,i], self.cosmics[:,:,i], self.noiseless[:,:,i], self.unmitigated[:,:,i] = p, c, n, u

			self.background = self.subcube.ccd.backgroundimage*self.ninstack
			self.noise = self.subcube.ccd.noiseimage*np.sqrt(self.ninstack)
			self.catalog = self.subcube.camera.catalog



		'''
        # set the cadence to this one
        self.camera.setCadence(self.cadence)
        np.save(os.path.join(self.camera.directory, 'observationdictionary.npy'), self.inputs)

        # reset the counter
        self.camera.counter = 0

        c = self.camera.ccd


        # expose this CCD, pulling out the images
        self.camera.ccd.expose(advancecounter=False, **self.inputs['expose'])

		self.photons, self.cosmics, self.noiseless, self.unmitigated

		# advance the counter by hand
        c.camera.advanceCounter()
		'''



	# is this necessary??????????
	def bin(self, nsubexposures=60, strategy=Central(n=10), plot=False):
		'''Bin together 2-second exposures, using some cosmic strategy.'''

		# make sure that we're starting with a 2-second exposure
		assert(self.cadence == self.subexposurecadence)

		# create an empty binned cube object that has the right size and shape
		binned = Cube(subject=self.subject, size=self.size, cadence=self.cadence*nsubexposures, n=(np.int(self.n/nsubexposures)) )

		# add an additional array to keep track of what the unmitigated lightcurves would look like
		binned.unmitigated = np.zeros_like(binned.photons)

		# loop over the x and y pixels
		self.speak('binning {0} cube by {1} subexposures into a {2} cube'.format(self.shape, nsubexposures, binned.shape))
		for x in np.arange(self.shape[0]):
			for y in np.arange(self.shape[1]):
				self.speak('   {x}, {y} out of ({size}, {size})'.format(x=x, y=y, size=self.size))
				timeseries = Timeseries1D(self, (x,y), nsubexposures=nsubexposures)
				strategy.calculate(timeseries)
				if plot:
					strategy.plot()
				binned.photons[x,y] = strategy.binned['flux']
				binned.cosmicinputs[x,y] = strategy.binned['naive'] - strategy.binned['nocosmics']
				binned.unmitigated[x,y] = strategy.binned['naive']
		return binned



	def display(self, cube=None, name='cube', limit=50):
		'''Use ds9 to display the image cube.'''
		self.speak('displaying (up to {0:.0f} exposures) of the pixel cube with ds9'.format(limit))
		if cube is None:
			cube = self.photons
		try:
			self.ds9
		except:
			self.ds9 = ds9(name)
		self.ds9.many(cube, limit=limit)

	@property
	def directory(self):
		'''Return path to a directory where this cube's data can be stored.'''
		dir = self.ccd.directory + 'cubes/'
		craftroom.utils.mkdir(dir)
		return dir

	@property
	def filename(self):
		'''Return a filename for saving/loading this cube.'''

		return self.directory + 'cube_{n:.0f}exp_at{cadence:.0f}s_{stacker}.npy'.format(n=self.n, cadence=self.cadence, stacker=self.stacker.name.replace(' ', ''))



	def save(self):
		'''Save this cube a 3D numpy array (as opposed to a series of FITS images).'''
		self.speak( "Saving cube to " + self.filename)
		np.save(self.filename, (self.photons, self.cosmics, self.noiseless, self.unmitigated, self.background, self.noise, self.catalog))

	def load(self, remake=False):
		'''Load this cube from a 3D numpy array, assuming one exists with the appropriate size, for this field, at this cadence, with the right number of exposures.'''
		self.speak("Trying to load simulated cube from " + self.filename)
		try:
			assert(remake==False)
			self.photons, self.cosmics, self.noiseless, self.unmitigated, self.background, self.noise, self.catalog = np.load(self.filename)
			self.speak( 'Loaded cube from ' + self.filename)
		except:
			self.speak('No saved cube was found; generating a new one.\n          (looked in {0})'.format( self.filename))
			self.simulate()
			self.save()

	def cubify(self, image):
		'''Slightly reshape an image, so it can be cast into operations on the whole cube.'''
		return image.reshape(self.xpixels, self.ypixels, 1)

	def median(self, which='photons'):
		'''The median image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'median'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.median(array, 2)
		return self.summaries[key+which]

	def mean(self, which='photons'):
		'''The median image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'mean'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.mean(array, 2)
		return self.summaries[key+which]


	def mad(self, which='photons'):
		'''The median of the absolute deviation image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'mad'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.median(np.abs(array - self.cubify(self.median(which))), 2)
		return self.summaries[key+which]

	def std(self, which='photons'):
		'''The standard deviation image.'''
		array = self.__dict__[which].astype(np.float64)
		key = 'std'
		try:
			self.summaries[key+which]
		except:
			self.summaries[key+which] = np.std(array, 2)
		return self.summaries[key+which]

	def sigma(self, which='photons', robust=True):
		if robust:
			return 1.4826*self.mad(which)
		else:
			return self.std(which)

	def oneWeirdTrickToTrimCosmics(self, threshold=4):
		'''To be used for simulating ground-based cosmics removal.'''
		shape = np.array(self.photons.shape)
		shape[-1] = 1
		bad = self.photons > (self.median() + threshold*self.sigma(robust=True)).reshape(shape)
		x, y, z = bad.nonzero()
		self.photons[x,y,z] = self.median()[x,y]
		self.speak('trimmed {0} cosmics from cube!'.format(np.sum(bad)))

	def master(self, which='photons'):
		'''The (calculated) master frame.'''
		return self.median(which)

	def nsigma(self, which='photons', robust=True):
		array = self.__dict__[which].astype(np.float64)
		return (array - self.cubify(self.median(which)))/self.cubify(self.sigma(which, robust=robust))

	def write(self, normalization='none'):
		'''Save all the images to FITS, inside the cube directory.'''

		# make a directory for the normalization used
		dir = self.directory + normalization + '/'
		craftroom.utils.mkdir(dir)

		if normalization == 'none':
			flux = self.photons
		if normalization == 'nsigma':
			flux = self.nsigma()

		# loop through the images in the cube
		for i in range(self.n):
			# pick some kind of normalization for the image
			image = flux[:,:,i]
			self.ccd.writeToFITS(image, dir + normalization + '_{0:05.0f}.fits'.format(i))


	def plot(self, normalization='none', ylim=None):

		# choose how to normalize the lightcurves for plotting
		if normalization.lower() == 'none':
			normalizationarray = 1.0
			ylabel = 'Photons'
		elif normalization.lower() == 'master':
			normalizationarray = self.master().reshape(self.xpixels, self.ypixels, 1)
			ylabel='Relative Flux'
		elif normalization.lower() == 'median':
			normalizationarray = self.median().reshape(self.xpixels, self.ypixels, 1)
			ylabel='Relative Flux'

		# create a relative light curve (dF/F, in most cases)
		photonsnormalized = self.photons/normalizationarray
		cosmicsnormalized = self.cosmics/normalizationarray
		noiselessnormalized = self.noiseless/normalizationarray

		something = cosmicsnormalized > 0
		cosmicsnormalized[something] += noiselessnormalized[something]

		# set up a logarithmic color scale (going between 0 and 1)
		def color(x):
			zero = np.min(np.log(self.master()*0.5))
			span = np.max(np.log(self.master())) - zero
			#normalized = (np.log(x) -  zero)/span
			normalized = (np.log(x) -  zero)/span
			return plt.matplotlib.cm.Blues_r(normalized)
			#return plt.matplotlib.cm.YlGn(normalized)

		# create a plot
		scale = 1.5
		plt.figure('{} | normalization={}'.format(self, normalization),
					figsize = (np.minimum(self.xpixels*scale,10),np.minimum(self.ypixels*scale, 10)), dpi=72)
		gs = plt.matplotlib.gridspec.GridSpec(self.xpixels,self.ypixels, wspace=0, hspace=0)

		# loop over pixels (in x and y directions)
		for i in range(self.ypixels):
			for j in range(self.xpixels):

				# set up axis sharing, so zooming on one plot zooms the other
				try:
					sharex, sharey = ax, ax
				except:
					sharex, sharey = None, None
				ax = plt.subplot(gs[-(i+1),j], sharex=sharex, sharey=sharey)

				# color the plot panel based on the pixel's intensity
				ax.patch.set_facecolor(color(self.median()[i,j]))


				ax.plot(photonsnormalized[i,j,:], color='black')
				ax.plot(noiselessnormalized[i,j,:], color='blue', alpha=0.5)
				ax.plot(cosmicsnormalized[i,j,:], color='red', alpha=0.8)

				if i == 0 and j == 0:
					plt.setp(ax.get_xticklabels(), rotation=90)
					ax.set_xlabel('Time')
					ax.set_ylabel(ylabel)
				else:
					plt.setp(ax.get_xticklabels(), visible=False)
					plt.setp(ax.get_yticklabels(), visible=False)


		if ylim is None:
			ax.set_ylim(np.min(photonsnormalized), np.max(photonsnormalized))
		else:
			ax.set_ylim(*ylim)
		plt.draw()

'''
Create some timeseries.
They can be either 1D (imagining a perfect, idealized photometer) or 3D
(like a timeseries of real 2D images, with stellar flux on multiple pixels).
'''

from imports import *

class Timeseries(Talker):
	'''Object to store a light curve, both unbinned and binned.'''

	def __init__(self, **kwargs):
		'''
		Initialize timeseries object, either as a toy model or from a pixel drawn from a simulated image cube.
		'''

		# set up the Talker
		Talker.__init__(self)

		self.scale = 1

		# create the timeseries, based on the remaining inputs (either 1D or 3D)
		self.create(**kwargs)

	@property
	def shape(self):
		'''The shape of arrays in this timeseries.'''
		return (self.nexposures, self.nsubexposures)

	def plot(self):
		'''*Very* simple plot of this timeseries.'''
		plt.figure('unbinned')
		plt.cla()
		try:
			x = self.x % self.tm.planet.period.value*60*60*24.0
		except AttributeError:
			x = self.x
		plt.plot(x.flatten(), self.flux.flatten(), linewidth=0, marker='.', alpha=0.3)
		plt.show();

	def __str__(self):
		'''Summary of this object.'''
		return "{nexposures} exposures, {nsubexposures} subexposures, cosmic rays {cosmicsamplitude:.2}X noise".format(**self.__dict__)


	def createTransit(self, cosmics=False):

		raise DeprecationWarning("Not sure if this works!")
		raise RuntimeError

		period = np.random.uniform()
		self.createSimple(noise=False, cosmics=False)
		p = craftroom.transit.Planet(period=period, t0=period/2.0)
		s = craftroom.transit.Star()
		i = craftroom.transit.Instrument()

		self.tm = craftroom.transit.TM(planet=p, star=s, instrument=i)
		self.tlc = craftroom.transit.TLC(self.x*self.exposurecadence/60.0/60.0/24.0, self.flux, self.subexposurenoise*np.ones_like(self.flux))
		self.tlc.linkModel(self.tm)
		self.flux = self.tm.model()
		self.addNoise()
		if cosmics:
			self.addCosmics()

class Timeseries1D(Timeseries):
	'''1D timeseries, set by a S/N and an amplitude of cosmic rays.'''
	def create(self, 	nexposures=1324,
						nsubexposures=900,
						subexposurecadence = 2.0,
						snr=100,
						amplitude=1.0,
						probabilityofcosmic = 0.001,
						):
		'''
		Initialize a 1D toy model timeseries, given a S/N and an amplitude and rate for cosmics.

		The parameters for a toy model timeseries are:
			nexposures=[how many binned exposures?]
			nsubexposures=[how many subexposures within each exposure?]
			subexposurecadence=[what's the cadence, in seconds, of the subexposures?]
			snr=[what's the S/N of the timeseries?]
			amplitude=[what's the amplitude of cosmic rays relative to binned exposure noise?],
			probabilityofcosmic=[what's the probability of a cosmic ray hit in a given sub-exposure?],
		'''

		# what's the total number of binned exposures (e.g., how many half-hour FFIs?)
		self.nexposures = nexposures

		# what's the number of subexposures per exposure?
		self.nsubexposures = nsubexposures

		# what's the cadence (in real seconds) of those subexposures?
		self.subexposurecadence = subexposurecadence

		# what's the S/N of the subexposure timeseries? (the flux will default to 1)
		self.snr = snr

		# what's the noise of a subexposure
		self.subexposurenoise = 1.0/self.snr

		# is this the right way to set these up???
		self.exposurenoise = self.subexposurenoise/np.sqrt(self.nsubexposures)

		# the cosmic ray amplitude is relative to the binned exposure noise
		self.cosmicsamplitude = amplitude

		# what's the rate of cosmic rays per binned exposure?
		self.cosmicspersubexposure = probabilityofcosmic
		self.cosmicsperexposure = self.cosmicspersubexposure*self.nsubexposures #1.0 - (1.0 - probabilityofhit)**self.nsubexposures
		self.exposurecadence = self.nsubexposures*self.subexposurecadence

		# then create a simple timeseries
		self.createSimple()


	def createSimple(self, cosmics=True, noise=True):
		'''Populate this timeseries with a simple toy model.'''

		# create an array to store flux values
		self.flux = np.ones(self.shape)

		# add noise to the timeseries
		if noise:
			self.addNoise()

		# add cosmic rays to the timeseries
		if cosmics:
			self.addCosmics()

		# define an axis of cadence numbers along the timeseries
		self.x = np.arange(0, self.nexposures, 1.0/self.nsubexposures).reshape(self.shape)

		# note that this is a toy model timeseries
		self.toy = True

	def addNoise(self):
		'''For toy model, include Gaussian noise.'''

		# make a noiseless timeseries
		self.noiseless = self.flux + 0.0

		# add a noise realization to each subexposure
		self.flux += np.random.normal(0,self.subexposurenoise,self.shape)

	def addCosmics(self):
		'''For toy model, include cosmic rays as random impulses with fixed amplitude.'''

		# a flux array with amplitude of the cosmic ray amplitude times the number of cosmics per subexposure
		self.cosmics = self.cosmicsamplitude*self.exposurenoise*self.nsubexposures*np.random.poisson(self.cosmicspersubexposure, self.shape)

		# add these cosmics into the timeseries
		self.flux += self.cosmics


class Timeseries1DTESS(Timeseries1D):
	'''1D cartoon timeseries, simulating a single TESS pixel.'''

	def __init__(self, 	nexposures=1324,
						nsubexposures=900,
						probabilityofcosmic=0.001,
						magnitude=10.0,
						subexposurecadence = 2.0,
						):
		'''
		Initialize a 1D toy model timeseries for a single TESS pixel timeseries,
		given a magnitude for the star, and the photons/cosmic ray.

		The parameters for a toy model timeseries are:
			nexposures=[how many binned exposures?]
			nsubexposures=[how many subexposures within each exposure?]
			probabilityofcosmic=[what's the probability of a cosmic ray hit in a given sub-exposure?],
			magnitude=[the brightness of the star sets the snr of the timeseries (and the amplitude)]
			)
		'''


		# what's the TESS magnitude of the star?
		self.magnitude = magnitude

		# what fraction of the stellar light is contained in that single pixel?
		self.containedflux = 0.3

		# how many stellar photons land on this pixel (per subexposure)?
		self.photonsfromstar = 1.7e6*subexposurecadence*73.0*10**(-0.4*self.magnitude)*self.containedflux

		# what's the noise from the sky, and from readout?
		skynoise = 10.0
		readoutnoise = 10.0

		# total number of stellar photons on the pixel in an exposure (in photoelectrons)
		signal = nsubexposures*self.photonsfromstar

		# calculate the noise (in photoelectrons)
		noise = np.sqrt(nsubexposures*(readoutnoise**2 + skynoise**2 + self.photonsfromstar))
		snr = signal/noise

		# calculate the cosmic ray impact, fractionally, compared to the binned exposure
		photonspercosmic=1000.0

		# this is assuming all cosmic rays are the same in their amplitude
		amplitude = photonspercosmic/noise

		Timeseries1D.__init__(self, nexposures=nexposures,
									nsubexposures=nsubexposures,
									subexposurecadence=subexposurecadence,
									snr=snr,
									amplitude=amplitude,
									probabilityofcosmic=probabilityofcosmic)

class Timeseries3D(Timeseries):

	def create(self, cube, pixel=(0,0)):
		#
		'''
		Initalize a timeseries from a particular pixel within an image cube.

				t = timeseries(cube=[an (pix)x(pix)x(time) image cube object], pixel=[tuple indicated which pixel to use, like '(4,3)'])
		'''
		self.createFromCube(cube, pixel)


	def createFromCube(self, cube, pixel):
		'''Populate this timeseries, using a pixel from a cube.'''

		self.nexposures = (np.int(cube.n/self.nsubexposures))
		x = pixel[0]
		y = pixel[1]
		# KLUDGE (possibly?) -- make sure that reshape is working the way we think it is
		self.flux = cube.photons[x,y,:self.nexposures*self.nsubexposures].reshape(self.shape)
		self.cosmics = cube.cosmics[x,y,:self.nexposures*self.nsubexposures].reshape(self.shape)
		self.subexposurenoise = cube.sigma()[x,y]
		self.exposurenoise = self.subexposurenoise/np.sqrt(self.nsubexposures)
		iscosmic = cube.cosmics.flatten() > 0
		self.cosmicsamplitude = np.mean(cube.cosmics.flatten()[iscosmic])/self.nsubexposures/self.exposurenoise
		self.cosmicspersubexposure = np.sum(iscosmic).astype(np.float)/len(cube.cosmics.flatten())
		self.cosmicsperexposure = self.cosmicspersubexposure*self.nsubexposures

		self.x = np.arange(0, self.nexposures, 1.0/self.nsubexposures).reshape(self.shape)
		self.toy = False
		self.cube = cube
		self.pixel = pixel

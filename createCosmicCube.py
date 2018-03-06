from Cube import *
from imports import *


"""
# create a cube, using these inputs
a = Cube(subarray=100, n=3, cadence=120)
a.load()
a.display(a.photons)


b = Cube(subarray=6, n=1, cadence=120)
b.simulate()
for n in ['master']:
    b.subcube.plot(normalization=n)
    plt.savefig('normaljitter_{}.pdf'.format(n))
"""

cosmicinputs['jitter']['amplifyinterexposurejitter'] = 1.0
normal = Cube(subarray=6, n=1, cadence=120, inputs=cosmicinputs)
normal.simulate()
for n in ['none', 'master']:
    normal.subcube.plot(normalization=n)
    plt.savefig('normaljitter_{}.pdf'.format(n))

cosmicinputs['jitter']['amplifyinterexposurejitter'] = 10.0
super = Cube(subarray=6, n=1, cadence=120, inputs=cosmicinputs)
super.simulate()
for n in ['none', 'master']:
    super.subcube.plot(normalization=n)
    plt.savefig('superjitter_{}.pdf'.format(n))

'''
a.save()

a.display(a.photons[:,:,0])
a.ds9.one(a.cosmics[:,:,0])
a.ds9.one(a.noiseless[:,:,0])
a.ds9.one((a.photons-a.noiseless)[:,:,0])

c = Cube(subarray=100, n=3, cadence=120)
c.load()




#plt.ion()
#c.subcube.plot(normalization='master')
#c.plot()
'''

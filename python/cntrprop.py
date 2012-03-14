from numpy import pi, zeros, arange, meshgrid, exp, sqrt 
from numpy import arctan as atan
L = 1.0
w0 = 1.0
sigma = 400
IL = 0.565
eta = 1.0
transits = 20

N = 256 # transverse grid slices
NZ = 20 # longitudinal grid slices
pumpoffset = [0.0, 0.0] # x then y
kin = 49.9
bkin = [0.0, 0.0] # backward k_in

noiseamp = 1e-7
seed = 12345678
probeamp = 2e-6
probewidth = 0.2
switchon = 60

relaxTransits = 5
transitsCompleted = 0
integrateSteps = (NZ-1)*(relaxTransits + transits)
zstep = L/(NZ-1)
xstep = 4.0*w0/N
ystep = 0.99*xstep
kxstep = 2*pi/(N*xstep)
kystep = 2*pi/(N*ystep)
lbynzeta = L / NZ * eta
stpp = 1j*lbynzeta
ramp = 1.0

f = zeros((N,N,NZ))
b = zeros((N,N,NZ))
f_tmp = zeros((N,N,NZ))
b_tmp = zeros((N,N,NZ))
ff = zeros((N,N))
bb = zeros((N,N))
expDz = zeros((N,N))
expDzByTwo = zeros((N,N))
noise = zeros((N,N))

x = arange(-N/2, N/2, 1)
y = arange(-N/2, N/2, 1)
i, j = meshgrid(x, y)

theta_n = zstep/(2*sigma) * ((i*kxstep)**2 + (j*kystep)**2)
expDz = exp(-1j * theta_n)
expDzByTwo = exp(-1j * theta_n/2.0)

def spot_size (z, z0, w0):
	return w0 * sqrt(1 + (z/z0)*(z/z0) )

def radius_curvature (z, z0):
	# This could be smarter, just adding epsilon to avoid nan's
	if (z == 0): 
		z += 1e-31
	return z * ( 1 + (z0/z)*(z0/z) )

def guoy_phase (z, z0):
	return atan(z/z0)

def rayleigh_range (w0, wavelambda):
	return pi*w0*w0/wavelambda;

def gaussian_beam (x, y, z, z0, w0, k):
	r = sqrt(x*x + y*y)
	w = spot_size(z,z0,w0)
	R = radius_curvature(z,z0)
	eta = guoy_phase(z,z0)
  # return w0/w * exp(- r*r/(w*w)) * exp(-I*k*z - I*k*r*r/(2*R) + I*eta);
  # The exp(ikz) term in the previous definition, causes extra phase accumulation
  # compared to the BPM. I need to sort this out for sure. Right now, leave it off.
	return w0/w * exp(- r*r/(w*w)) * exp(- 1j*k*r*r/(2*R) + 1j*eta)



# do something to add noise

for intstep in range(integrateSteps):
	print "integration step %d of %d" % (intstep, integrateSteps)
	f[:,:,0] = sqrt(IL) * ramp * gaussian_beam(x - pumpoffset[0], y - pumpoffset[1], 0, L*sigma/2, w0, L*sigma/(w0*w0)) * exp(1j*(bkin[0]*x + bkin[1]*y)) + sqrt(IL)*ramp*probeamp*exp(- (x*x + y*y) / ((w0*probewidth)*(w0*probewidth)))*exp(1j*kin*x) # + ramp*noiseamp*noise(i,j) #TODO restore noise term
	b[:,:,NZ-1] = sqrt(IL) * ramp * gaussian_beam(x+pumpoffset[0], y+pumpoffset[1], 0, L*sigma/2, w0, L*sigma/(w0*w0))



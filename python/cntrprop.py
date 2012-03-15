from numpy import pi, zeros, arange, meshgrid, exp, sqrt, conjugate 
from numpy.fft import fftshift, fft2, ifft2
from numpy import arctan as atan
from sys import stdout
L = 1.0
w0 = 1.0
sigma = 400
IL = 0.565
eta = 1.0
transits = 1

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

f = zeros((N,N,NZ)).astype(complex) # the fwd field at all points
b = zeros((N,N,NZ)).astype(complex) # bwd field at all points
f_tmp = zeros((N,N,NZ)).astype(complex) # tmp version of fwd
b_tmp = zeros((N,N,NZ)).astype(complex) # tmp version of bwd
probe = zeros((N,N)).astype(complex) # probe field
ff = zeros((N,N)).astype(complex) # use to calculate new slice of f
bb = zeros((N,N)).astype(complex) # use to calculate new slice of b
expDz = zeros((N,N)).astype(complex) # full-step free-space propagator
expDzByTwo = zeros((N,N)).astype(complex) # half-step FS propagator
noise = zeros((N,N)) # a noise term

x = arange(-N/2, N/2, 1)
y = arange(-N/2, N/2, 1)
i, j = meshgrid(x, y)
x = i*xstep # set master coordinate grid
y = j*ystep

theta_n = zstep/(2*sigma) * ((i*kxstep)**2 + (j*kystep)**2)

expDz = fftshift(exp(-1j * theta_n)) # fftshift to apply in FT space
expDzByTwo = fftshift(exp(-1j * theta_n/2.0))

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

probe = sqrt(IL)*ramp*probeamp*exp(- (x*x + y*y) / ((w0*probewidth)*(w0*probewidth)))*exp(1j*kin*x)

# do something to add noise

for intstep in range(integrateSteps):
    print "integration step %d of %d" % (intstep, integrateSteps)
    f[:,:,0] = sqrt(IL) * ramp * gaussian_beam(x - pumpoffset[0], y - pumpoffset[1], 0, L*sigma/2, w0, L*sigma/(w0*w0)) * exp(1j*(bkin[0]*x + bkin[1]*y))# + ramp*noiseamp*noise(i,j) #TODO restore noise term
    if intstep > switchon: 
        f[:,:,0] += probe
    b[:,:,NZ-1] = sqrt(IL) * ramp * gaussian_beam(x+pumpoffset[0], y+pumpoffset[1], 0, L*sigma/2, w0, L*sigma/(w0*w0))

    # step forward at the first point
    # first find the intensities of each field:
    fi = f[:,:,0]*conjugate(f[:,:,0])
    bi = b[:,:,0]*conjugate(b[:,:,0])
    # then apply the nonlinear operator to first slice:
    ff = f[:,:,0] * exp(stpp*(fi+2.0*bi))

    # step backward at the last points:
    # intensities at the far end
    fi = f[:,:,-1]*conjugate(f[:,:,-1])
    bi = b[:,:,-1]*conjugate(b[:,:,-1])
    # nonlinear operator on last slice:
    bb = b[:,:,-1] * exp(stpp*(2.0*fi+bi))

    # Fourier transform the end slices
    ff = fft2(ff)
    bb = fft2(bb)
    print "f"

    # Free-space propagate half-step for end slices
    ff = ff*expDzByTwo
    bb = bb*expDzByTwo
    print "."

    # inverse Fourier transform end slices
    ff = ifft2(ff)
    bb = ifft2(bb)
    print "b"

    # set field values one step in from each end (use tmp array)
    f_tmp[:,:,1] = ff
    b_tmp[:,:,-2] = bb

    # LOOP: for all slices
    for zindex in range(1,NZ-1):
        print "."
        
        # nonlinear step for each slice
        fi = f[:,:,zindex]*conjugate(f[:,:,zindex]) 
        bi = b[:,:,zindex]*conjugate(b[:,:,zindex]) 
        ff = f[:,:,zindex] * exp(stpp*(fi + 2.0*bi))
        bb = b[:,:,zindex] * exp(stpp*(2.0*fi + bi))
        
        # Forward FFT
        ff = fft2(ff)
        bb = fft2(bb)

        # TODO: save if it is the last slice

        # Free-space propagate a full step
        ff *= expDz
        bb *= expDz

        # Backward FFT
        ff = ifft2(ff)
        bb = ifft2(bb)
        
        # set field values in next (and prev) slice for fwd (bwd) fields
        f_tmp[:,:,zindex+1] = ff
        b_tmp[:,:,zindex-1] = bb
        
    # shift back from working version.
    f = f_tmp
    b = b_tmp

    if (intstep > 1) and ((intstep + 1) % (NZ-1) == 0):
        transitsCompleted = (intstep+1)/(NZ-1)
        print "transits: " + str(transitsCompleted)

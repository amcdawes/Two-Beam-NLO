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
integratesteps = (NZ-1)*(relaxTransits + transits)
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

x = np.arange(-N/2, N/2, 1)
y = np.arange(-N/2, N/2, 1)
i, j = np.meshgrid(x, y)

theta_n = zstep/(2*sigma) * ((i*kxstep)**2 + (j*kystep)**2)
expDz = exp(-1j * theta_n)
expDzByTwo = exp(-1j * theta_n/2.0)

# do something to add noise

for intstep in range(integrateSteps):
	print "integration step %d of %d" % (intstep, integrationSteps)
	f = sqrt(IL) * ramp * gaussian_beam(x - pumpoffset[0], y - pumpoffset[1], 0, L*sigma/2, w0, L*sigma/(w0*w0) * exp(1j*(bkin[0]*x + bkin[1]*y)) + sqrt(IL)*ramp*


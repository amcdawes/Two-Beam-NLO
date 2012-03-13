// Gaussian beam parameters, put this in a header later:

double spot_size (double z, double z0, double w0)
{
  return w0 * sqrt(1 + (z/z0)*(z/z0) );
}

double radius_curvature (double z, double z0)
{
  if (z == 0) // This could be smarter, just adding epsilon to avoid nan's
  {
    z += 1e-31;
  }
  return z * ( 1 + (z0/z)*(z0/z) );
}

double guoy_phase (double z, double z0)
{
  return atan(z/z0);
}

double rayleigh_range (double w0, double lambda)
{
  return M_PI*w0*w0/lambda;
}

Complex gaussian_beam (double x, double y, double z, double z0, double w0, double k) 
{
  Complex I(0,1);
  double r = sqrt(x*x + y*y);
  double w = spot_size(z,z0,w0);
  double R = radius_curvature(z,z0);
  double eta = guoy_phase(z,z0);
  // return w0/w * exp(- r*r/(w*w)) * exp(-I*k*z - I*k*r*r/(2*R) + I*eta);
  // The exp(ikz) term in the previous definition, causes extra phase accumulation
  // compared to the BPM. I need to sort this out for sure. Right now, leave it off.
  return w0/w * exp(- r*r/(w*w)) * exp(- I*k*r*r/(2*R) + I*eta);
}

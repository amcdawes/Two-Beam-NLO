#include <iostream>
#include <cmath>
#include <pngwriter.h>

#include "Array.h"
#include "fftw++.h"
#include "gaussian.h"

#define DEBUG 0
#define DUMP 0
#define DATA 0

using namespace std;
using Array::array3;
using Array::array2;

void array_to_png (char* filename, int size, array2<Complex> &field, int frame)
{
  char text[30];
  int i,j;
  double pixel;
  pngwriter png(size,size,0,filename);
  // Peak-finding for image scaling:
  double max = 1e-32;
  for (i = 0; i < size; i++)
  {
    for (j = 0; j < size; j++)
    {
      if (real(field(i,j)*conj(field(i,j))) > max) max = real(field(i,j)*conj(field(i,j)));
    }
  }
  
  // Image writeout.
  for (i = 0; i < size; i++)
  {
    for (j = 0; j < size; j++)
    {
      pixel = real(field(i,j)*conj(field(i,j)));
      png.plot(i,j,pixel/max,pixel/max,0.0);
    }
  }

  sprintf(text, "%i", frame);
  png.plot_text("/xtmp/dawes/share/pngwriter/fonts/FreeMonoBold.ttf",10,10,10,0.0,text,1.0,1.0,1.0);
  png.close();
}

void threed_array_to_png (char* filename, int size, int length, array3<Complex> &field, int frame)
{
  char text[30];
  int i,j;
  double pixel;
  pngwriter png(size,size,0,filename);
  // Peak-finding for image scaling:
  double max = 1e-32;
  for (i = 0; i < size; i++)
  {
    for (j = 0; j < size; j++)
    {
      if (real(field(i,j,length-1)*conj(field(i,j,length-1))) > max) max = real(field(i,j,length-1)*conj(field(i,j,length-1)));
    }
  }
  
  //cout << max1 << "\t" << max2 << endl;
  
  // Image writeout.
  for (i = 0; i < size; i++)
  {
    for (j = 0; j < size; j++)
    {
      pixel = real(field(i,j,length-1)*conj(field(i,j,length-1)));
      png.plot(i,j,pixel/max,pixel/max,0.0);
    }
  }
  sprintf(text, "%i", frame);
  png.plot_text("/xtmp/dawes/share/pngwriter/fonts/FreeMonoBold.ttf",10,10,10,0.0,text,1.0,1.0,1.0);
  png.close();
}

void save_intensity(char* filename, int NX, int NY, int slice, array3<Complex> &field)
{
	int i, j;
	ofstream dumpfile(filename);
	for (i = 0; i < NX; i++)
	{
		for (j = 0; j < NY; j++)
		{
      dumpfile << real(field(i,j,slice)*conj(field(i,j,slice))) << "\t";
      //dumpfile << imag(conj(field(i,j,slice))) << "\t";
		}
		dumpfile << endl;
	}
}

int
main (void)
{
  Complex I(0,1);
  int i, j, k, n, m, intstep, zindex;
  double x, y;

  // PHYSICAL PARAMETERS
  // I use L and w0 later even though they are both 1. This is for clarity 
  // and readability.
  double L = 1.0; // This MUST be 1.0, the simulation is scaled as such
  double w0 = 1.0; // This also must be 1.0!
  double sigma = 400.0; // This is a knob for changing the scale.
  double IL = 0.565; // Unitless Intensity
  double eta = 1.0;
  int transits = 20;

  
  // User input:
  cout << "Enter sigma (1-1000): ";
  cin >> sigma;
  cout << "Enter transits (1-100): ";
  cin >> transits;
  cout << "Enter IL (0.4 - 5.0): ";
  cin >> IL;
  cout << "sigma: " << sigma << " transits: " << transits << " IL: " << IL << endl;
  
  // NUMERICAL PARAMETERS (GRID)
  int N = 256; // Transverse grid points
  int NZ = 20; // Longitudinal grid points
  int relaxTransits = 5;
  int transitsCompleted = 0;
  int integrateSteps = (NZ-1)*(relaxTransits + transits);
  char imagefile[30];
  char outfieldfile[30];
  double zstep = L/(NZ-1);
  double xstep = 6.0*w0/N;
  double ystep = xstep;
  double kxstep = 2*M_PI/(N*xstep);
  double kystep = 2*M_PI/(N*ystep);
  double lbynzeta = L / NZ * eta;
  Complex stpp = I * lbynzeta;
  cout << "stpp: " << stpp << endl;
  double theta_n;
  double kin = 0.0;
  double ramp = 1.0; // use for slowly turning on fields
  double noiseamp = 0.00001;
  Complex fi(0,0), bi(0,0);
  
  size_t align = sizeof(Complex);
  
  array3<Complex> f(N,N,NZ,align);
  array3<Complex> b(N,N,NZ,align);
  array3<Complex> f_tmp(N,N,NZ,align);
  array3<Complex> b_tmp(N,N,NZ,align);
  array2<Complex> ff(N,N,align);
  array2<Complex> bb(N,N,align);
  array2<Complex> farfield(N,N,align);
  array2<Complex> expDz(N,N,align);
  array2<Complex> expDzByTwo(N,N,align);
  array2<Complex> noise(N,N,align);
  
  srandom(18765432);
  
  fft2d ffForward2(-1,ff);
  fft2d ffBackward2(1,ff);
  fft2d bbForward2(-1,bb);
  fft2d bbBackward2(1,bb);
  
  if (DEBUG) cerr << "Initialized" << endl;
  
  // The free-space propagator is ugly without the nice fftshift, but whatever.
  // now it also includes spatial filtering (theta_n < pi/2)
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      if (i < N/2 && j < N/2)
      {
      	theta_n = (zstep/(2*sigma)*(((i)*kxstep*(i)*kxstep) + ((j)*kystep*(j)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDz(i,j) = 0;
      	} else { expDz(i,j) = exp(-I*theta_n); }
      }
      else if (i < N/2 && j >= N/2)
      {
        theta_n = (zstep/(2*sigma)*(((i)*kxstep*(i)*kxstep) + ((j-N)*kystep*(j-N)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDz(i,j) = 0;
      	} else { expDz(i,j) = exp(-I*theta_n); }
      }
      else if (i >= N/2 && j < N/2)
      {
        theta_n = (zstep/(2*sigma)*(((i-N)*kxstep*(i-N)*kxstep) + ((j)*kystep*(j)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDz(i,j) = 0;
      	} else { expDz(i,j) = exp(-I*theta_n); }
      }
      else if (i >= N/2 && j >= N/2)
      {
        theta_n = (zstep/(2*sigma)*(((i-N)*kxstep*(i-N)*kxstep) + ((j-N)*kystep*(j-N)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDz(i,j) = 0;
      	} else { expDz(i,j) = exp(-I*theta_n); }
      } 
      else
      {
        cerr << "ERROR: could not properly determine the propagator" << endl;
      }
    }
  }
  
  // Now calculate expDzByTwo:
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      if (i < N/2 && j < N/2)
      {
      	theta_n = (zstep/(4*sigma)*(((i)*kxstep*(i)*kxstep) + ((j)*kystep*(j)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDzByTwo(i,j) = 0;
      	} else { expDzByTwo(i,j) = exp(-I*theta_n); }
      }
      else if (i < N/2 && j >= N/2)
      {
        theta_n = (zstep/(4*sigma)*(((i)*kxstep*(i)*kxstep) + ((j-N)*kystep*(j-N)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDzByTwo(i,j) = 0;
      	} else { expDzByTwo(i,j) = exp(-I*theta_n); }
      }
      else if (i >= N/2 && j < N/2)
      {
        theta_n = (zstep/(4*sigma)*(((i-N)*kxstep*(i-N)*kxstep) + ((j)*kystep*(j)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDzByTwo(i,j) = 0;
      	} else { expDzByTwo(i,j) = exp(-I*theta_n); }
      }
      else if (i >= N/2 && j >= N/2)
      {
        theta_n = (zstep/(4*sigma)*(((i-N)*kxstep*(i-N)*kxstep) + ((j-N)*kystep*(j-N)*kystep)));
      	if (theta_n > M_PI/2)
      	{
      		expDzByTwo(i,j) = 0;
      	} else { expDzByTwo(i,j) = exp(-I*theta_n); }
      } 
      else
      {
        cerr << "ERROR: could not properly determine the expDzByTwoByTwo propagator" << endl;
      }
    }
  }
  

  if (DUMP)
  {
  	ofstream dumpfile("dumpfile.dat");
  	for (i = 0; i < N; i++)
  	{
  		for (j = 0; j < N; j++)
  		{
  			dumpfile << real(expDz(i,j)) << "\t";
  		}
  		dumpfile << endl;
  	}
  }
  
  for (i = 0; i < N; i++) 
  {
    for (j = 0; j < N; j++)
    {
      noise(i,j) = double(random()) / double(RAND_MAX); // Set noise term.
    }
  }
  
  // Propagate to L
  if (DEBUG) cerr << "Beginning integration" << endl;
  for (intstep = 0; intstep < integrateSteps; intstep++)
  {
    cerr << "Starting step: " << intstep+1 << " out of " << integrateSteps << endl;

    // Set front and back of medium with ramp functions
    ramp = 1 - exp(- intstep*zstep );
    // cout << "Ramp = " << ramp << endl;
    for (i = 0; i < N; i++) 
    {
      for (j = 0; j < N; j++)
      {
        x = (i-N/2)*xstep;
        y = (j-N/2)*ystep;
        f(i,j,0) = sqrt(IL) * ramp * gaussian_beam(x, y, 0, L*sigma/2, w0, L*sigma/(w0*w0)) * exp(I*kin*(x+y)) + noiseamp*(1-ramp)*noise(i,j);
        //b_tmp(i,j,NZ-1) = sqrt(IL) * ramp * gaussian_beam(x, y, 0, L*sigma/2, w0, L*sigma/(w0*w0)) + noiseamp*noise(i,j);
		// Matlab is only adding noise to f right now.
        b(i,j,NZ-1) = sqrt(IL) * ramp * gaussian_beam(x, y, 0, L*sigma/2, w0, L*sigma/(w0*w0));      
      }
    }
    
    // first and last points individually 
    // nonlinear part
    for (i = 0; i < N; i++) 
    {
      for (j = 0; j < N; j++)
      {
        fi = f(i,j,0)*conj(f(i,j,0));
        bi = b(i,j,0)*conj(b(i,j,0));
        ff(i,j) = f(i,j,0)*exp(stpp*(fi+2.0*bi)); // forward at first point
        fi = f(i,j,NZ-1)*conj(f(i,j,NZ-1));
        bi = b(i,j,NZ-1)*conj(b(i,j,NZ-1));
        bb(i,j) = b(i,j,NZ-1)*exp(stpp*(2.0*fi+bi)); // backward at last point
      }
    }
    
    //cout << ff << endl;
    
    // Fourier transform
    ffForward2.fft(ff);
    bbForward2.fft(bb);
    
    if (DEBUG) cerr << "Forward FFT complete" << endl;
    
    // free-space propagate half the distance
    for (i = 0; i < N; i++) 
    {
      for (j = 0; j < N; j++)
      {
        ff(i,j) = ff(i,j) * expDzByTwo(i,j);
        bb(i,j) = ff(i,j) * expDzByTwo(i,j);
      }
    }
    
    // inverse Fourier transform
    ffBackward2.fftNormalized(ff);
    bbBackward2.fftNormalized(bb);
    
    if (DEBUG) cerr << "Backward FFT complete" << endl;
    
    //cout << ff << endl;
    
    // set field values one step in from each end
    for (i = 0; i < N; i++) 
    {
      for (j = 0; j < N; j++)
      {
        f_tmp(i,j,1) = ff(i,j); 
        b_tmp(i,j,NZ-2) = bb(i,j);
      }
    }
    
    if (DEBUG) cerr << "Integrating along Z" << endl;
    
    for (zindex = 1; zindex < NZ-1; zindex++) // integrate along z axis
    { 
      // This is the loop to parallelize. Each iteration can be individually
      // processed: input f(N,N) and b(N,N) at that slice and returning ff(N,N) 
      // and bb(N,N) which are then assigned to f_tmp and b_tmp in +1 and -1 pos.
      cerr << ".";
      // calculate intensity at these slices:
      // nonlinear propagate each transverse point
      for (i = 0; i < N; i++) 
      {
        for (j = 0; j < N; j++)
        {
          fi = f(i,j,zindex)*conj(f(i,j,zindex));
          bi = b(i,j,zindex)*conj(b(i,j,zindex));
          ff(i,j) = f(i,j,zindex)*exp(stpp*(fi + 2.0*bi));
          bb(i,j) = b(i,j,zindex)*exp(stpp*(2.0*fi + bi));
        }
      }
      
      // Fourier transform
      ffForward2.fft(ff);
      bbForward2.fft(bb);

      // if this is the last slice, save the fourier transform since we just calculated it.
      if (zindex == NZ-2 && ((intstep+1) % (NZ-1) == 0))
      {
        cout << "saving farfield" << endl;
        for (i = 0; i < N; i++) 
        {
          for (j = 0; j < N; j++)
          {
            x = (i-N/2)*xstep;
            y = (j-N/2)*ystep;
            if ((n = i - N/2) < 0) n+=N;
            if ((m = j - N/2) < 0) m+=N;
            farfield(i,j) = ff(n,m)*(1 - exp(-(x*x+y*y)/(0.55*0.55)));
          }
        }       
        sprintf(imagefile, "%03i_farfield.png", transitsCompleted);
        array_to_png(imagefile, N, farfield, (intstep+1)/(NZ-1));
      }

      // free-space propagate over the full step
      for (i = 0; i < N; i++) 
      {
        for (j = 0; j < N; j++)
        {
          ff(i,j) *= expDz(i,j);
          bb(i,j) *= expDz(i,j);
        }
      }

      // inverse Fourier transform
      ffBackward2.fftNormalized(ff);
      bbBackward2.fftNormalized(bb);
      
      // set field values in next (previous) z-plane for fwd (bwd) field
      for (i = 0; i < N; i++) 
      {
        for (j = 0; j < N; j++)
        {
          f_tmp(i,j,zindex+1) = ff(i,j); 
          b_tmp(i,j,zindex-1) = bb(i,j);
        }
      }
      
    }

    // shift fields back from working version
    for (i = 0; i < N; i++) 
    {
      for (j = 0; j < N; j++)
      {
        for (k = 0; k < NZ; k++)
        {
          f(i,j,k) = f_tmp(i,j,k);
          b(i,j,k) = b_tmp(i,j,k);
        }
      }
    }
    cerr << endl; // end the line of status dots
    if (intstep > 1 && (intstep+1) % (NZ-1) == 0) 
    {
         transitsCompleted = (intstep+1) / (NZ-1);
         cout << "transits: " << transitsCompleted << endl;
         sprintf(imagefile, "%03i_frame.png", transitsCompleted);
         sprintf(outfieldfile, "%03i_out.dat", transitsCompleted);
         //sprintf(infieldfile, "%03i_in.dat", transitsCompleted);         
         if (DATA) save_intensity(outfieldfile, N, N, NZ-1, f);
         //save_intensity(infieldfile, N, N, 0, f);
         threed_array_to_png(imagefile, N, NZ, f, transitsCompleted);
    }
  } // end of BPM
  
  
  // calculate far-field intensity pattern from the final forward slice.
  for (i = 0; i < N; i++) 
  {
    for (j = 0; j < N; j++)
    {
      ff(i,j) = f(i,j,NZ-1);
    }
  }
  
  ffForward2.fft(ff);
  for (i = 0; i < N; i++) 
  {
    for (j = 0; j < N; j++)
    {
      x = (i-N/2)*xstep;
      y = (j-N/2)*ystep;
      if ((n = i - N/2) < 0) n+=N;
      if ((m = j - N/2) < 0) m+=N;
      farfield(i,j) = ff(n,m)*(1 - exp(-(x*x+y*y)/(0.55*0.55)));
    }
  }
  
  array_to_png("farfield.png", N, farfield,0);
  
  ofstream outfile("datafile.dat");
  
  for (i = 0; i < N; i++) {
    outfile << real(f(i,N/2,NZ-1)*conj(f(i,N/2,NZ-1)))
    << "\t" << real(b(i,N/2,0)*conj(b(i,N/2,0)))
    << endl;
  }
  
  return 0;
}

#include <iostream>
#include <cmath>
#include <pngwriter.h>

#include "Array.h"
#include "fftw++.h"
#include "gaussian.h"
#include <omp.h>

#define DEBUG 0
#define NFDATA 0
#define FFDATA 1

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
  png.plot_text("/usr/X11/lib/X11/fonts/TTF/VeraMono.ttf",10,10,10,0.0,text,1.0,1.0,1.0);
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
  png.plot_text("/usr/X11/lib/X11/fonts/TTF/VeraMono.ttf",10,10,10,0.0,text,1.0,1.0,1.0);
  png.close();
}

void save_intensity(char* filename, int NX, int NY, int slice, array3<Complex> &field)
{
	int i, j;
	ofstream dumpfile(filename);
	for (j = 0; j < NY; j++)
	{
		for (i = 0; i < NX; i++)
		{
      dumpfile << real(field(i,j,slice)*conj(field(i,j,slice))) << "\t";
      //dumpfile << imag(conj(field(i,j,slice))) << "\t";
		}
		dumpfile << endl;
	}
}

void save_twod_complex_array(char* filename1, char* filename2, int NX, int NY, array2<Complex> &field)
{
	int i, j;
	ofstream realfile(filename1);
	ofstream imagfile(filename2);
	for (j = 0; j < NY; j++)
	{
		for (i = 0; i < NX; i++)
		{
      		realfile << real(field(i,j)) << "\t";
      		imagfile << imag(field(i,j)) << "\t";
		}
		realfile << endl;
		imagfile << endl;
	}
}

void save_twod_real_array(char* filename1, int NX, int NY, array2<Complex> &field)
{
	int i, j;
	ofstream realfile(filename1);
	for (j = 0; j < NY; j++)
	{
		for (i = 0; i < NX; i++)
		{
      		realfile << real(field(i,j)*conj(field(i,j))) << "\t";
		}
		realfile << endl;
	}
}

void save_twod_real_array_cropped(char* filename1, int NX, int NY, int start, int stop, array2<Complex> &field)
{
	int i, j;
	ofstream realfile(filename1);
	for (j = start; j < stop+1; j++)
	{
		for (i = start; i < stop+1; i++)
		{
      		realfile << real(field(i,j)*conj(field(i,j))) << "\t";
		}
		realfile << endl;
	}
}

int
main (void)
{
  Complex I(0,1);
  int i, j, k, n, m, intstep, zindex;
  double x, y;
  char imagefile[30];
  char outfieldfile[30];
  char outintfile[30];
  char noisefiler[30];
  char noisefilei[30];

  // PHYSICAL PARAMETERS
  // I use L and w0 later even though they are both 1. This is for clarity 
  // and readability.
  double L = 1.0; // This MUST be 1.0, the simulation is scaled as such!
  double w0 = 1.0; // This also must be 1.0!
  double sigma = 400.0; // This is a knob for changing the scale.
  double IL = 0.565; // Unitless Intensity
  double eta = 1.0;
  int transits = 20;
  
  // NUMERICAL PARAMETERS (grid size)
  int N = 256; // Transverse grid points
  int NZ = 20; // Longitudinal grid points
  double pumpoffsetx = 0.0;
  double pumpoffsety = 0.0;
  double kin = 49.0;
  double bkinx = 0.0;
  double bkiny = 0.0;
  double noiseamp = 0.0000001;
  int seed = 12345678;
  double probeamp = 0.000002;
  double probewidth = 0.2;
  int switchon = 60;
    
  // User input:
  cout << "Enter sigma (1-1000): ";
  cin >> sigma;
  cout << "Enter transits (1-100): ";
  cin >> transits;
  cout << "Enter IL (0.4 - 5.0): ";
  cin >> IL;
  cout << "Enter NZ (10-50): ";
  cin >> NZ;
  cout << "Enter pumpoffsetx (0.0): ";
  cin >> pumpoffsetx;
  cout << "Enter pumpoffsety (0.0): ";
  cin >> pumpoffsety;
  cout << "Enter kin (~12-49): ";
  cin >> kin;
  cout << "Enter backward k in x (0.0): ";
  cin >> bkinx;
  cout << "Enter backward k in y (0.0): ";
  cin >> bkiny;
  cout << "Enter noiseamp (0.0000001): ";
  cin >> noiseamp;
  cout << "Enter seed (12345678): ";
  cin >> seed;
  cout << "Enter probeamp (0.000002): ";
  cin >> probeamp;
  cout << "Enter probewidth (0.2): ";
  cin >> probewidth;
  cout << "Enter Tr for switch turn-on (60): ";
  cin >> switchon;
  
  ofstream logfile("c_log.dat"); // somewhat for debugging, but also to log the parameters
  
  logfile << "sigma: " << sigma 
  		  << " transits: " << transits 
  		  << " IL: " << IL 
  		  << " NZ: " << NZ 
  		  << " kin: " << kin
  		  << " bkinx: " << bkinx
  		  << " bkiny: " << bkiny
  		  << " noiseamp: " << noiseamp
  		  << " seed: " << seed
  		  << " probeamp: " << probeamp
  		  << " probewidth: " << probewidth
  		  << " switch-on at: " << switchon
  		  << endl;
 
 
  // NUMERICAL PARAMETERS (calculated)
  int relaxTransits = 5;
  int transitsCompleted = 0;
  int integrateSteps = (NZ-1)*(relaxTransits + transits);
  double zstep = L/(NZ-1);
  double xstep = 4.0*w0/N;
  double ystep = 0.99*xstep;
  double kxstep = 2*M_PI/(N*xstep);
  double kystep = 2*M_PI/(N*ystep);
  double lbynzeta = L / NZ * eta;
  double randoma = 0;
  double randomb = 0;
  Complex stpp = I * lbynzeta;
  cout << "stpp: " << stpp << endl;
  double theta_n;
  double ramp = 1.0; // use for slowly turning on fields

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
  
  srandom(seed);
  
  fft2d ffForward2(-1,ff);
  fft2d ffBackward2(1,ff);
  fft2d bbForward2(-1,bb);
  fft2d bbBackward2(1,bb);
  
  if (DEBUG) cerr << "Initialized" << endl;
  
  // The free-space propagator is ugly without the nice fftshift, but whatever.
  // now it also includes spatial filtering (theta_n < pi/2)
  //  AMCD March 5, changed to theta_n < pi to see if hexagons form at sigma = 35
  //  AMCD March 10, changed back to theta_n < pi/2 
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
        cerr << "ERROR: could not properly determine the expDzByTwo propagator" << endl;
      }
    }
  }
  

 
  // set the noise, dynamic noise in FWD direction
  for (i = 0; i < N; i++) 
  {
    for (j = 0; j < N; j++)
    {
      randoma = double(random()) / double(RAND_MAX); // first random number
      randomb = double(random()) / double(RAND_MAX); // second random number
      noise(i,j) = sqrt(-2*log(randoma)) * (cos(2*M_PI*randomb) + I * sin(2*M_PI*randomb)); // Box-Mueller transform
      //noise(i,j) = sqrt(-2*log(randoma)) * (cos(2*M_PI*randomb)); // Box-Mueller transform, real only
    }
  }


  // sprintf(noisefiler, "realnoise.dat");
  // sprintf(noisefilei, "imagnoise.dat");
  // save_twod_complex_array(noisefiler, noisefilei, N, N, noise); // save the noise data to check
  
  // Propagate to L
  if (DEBUG) cerr << "Beginning integration" << endl;
  for (intstep = 0; intstep < integrateSteps; intstep++)
  {
    cerr << "Starting step: " << intstep+1 << " out of " << integrateSteps << endl;

    // Set front and back of medium with ramp functions
    ramp = 1 - exp(- intstep*zstep/5 );
    // cout << "Ramp = " << ramp << endl;

    if ( (transitsCompleted/switchon)%3 == 2) // apply the probe beam
    {
      for (i = 0; i < N; i++) 
      {
        for (j = 0; j < N; j++)
        {
          x = (i-N/2)*xstep;
          y = (j-N/2)*ystep;
          f(i,j,0) = sqrt(IL) * ramp * gaussian_beam(x-pumpoffsetx, y-pumpoffsety, 0, L*sigma/2, w0, L*sigma/(w0*w0)) * exp(I*(bkinx*x + bkiny*y)) + sqrt(IL)*ramp*probeamp*exp(- (x*x + y*y) / ((w0*probewidth)*(w0*probewidth)))*exp(I*kin*x) + ramp*noiseamp*noise(i,j);
          b(i,j,NZ-1) = sqrt(IL) * ramp * gaussian_beam(x+pumpoffsetx, y+pumpoffsety, 0, L*sigma/2, w0, L*sigma/(w0*w0));     
        }
      }
      
    } else { // inject only pump+noise 
			for (i = 0; i < N; i++) 
      {
        for (j = 0; j < N; j++)
        {
          x = (i-N/2)*xstep;
          y = (j-N/2)*ystep;
          f(i,j,0) = sqrt(IL) * ramp * gaussian_beam(x-pumpoffsetx, y-pumpoffsety, 0, L*sigma/2, w0, L*sigma/(w0*w0)) * exp(I*(bkinx*x + bkiny*y)) + ramp*noiseamp*noise(i,j);
          b(i,j,NZ-1) = sqrt(IL) * ramp * gaussian_beam(x+pumpoffsetx, y+pumpoffsety, 0, L*sigma/2, w0, L*sigma/(w0*w0));      
        }
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
    
	#pragma omp parallel for private(fi,bi,ff,bb)
    for (zindex = 1; zindex < NZ-1; zindex++) // integrate along z axis
    { 
      // This is the loop that is parallelized. Each iteration can be individually
      // processed: input f(N,N) and b(N,N) at that slice and returning ff(N,N) 
      // and bb(N,N) which are then assigned to f_tmp and b_tmp in +1 and -1 pos.
      // Each thread will handle a few of the steps in z, and at the end, the
      // entire f_tmp and b_tmp array will be constructed. Once it is ready,
      // it is just copied over f and b.

      

      // calculate intensity at these slices:
      // then compute the nonlinear propagate of each transverse point
	  int ii,jj,kk; // new local (private) iterators for parallel code
	  for (ii = 0; ii < N; ii++) 
      {
        for (jj = 0; jj < N; jj++)
        {
          fi = f(ii,jj,zindex)*conj(f(ii,jj,zindex));
          bi = b(ii,jj,zindex)*conj(b(ii,jj,zindex));
          ff(ii,jj) = f(ii,jj,zindex)*exp(stpp*(fi + 2.0*bi));
          bb(ii,jj) = b(ii,jj,zindex)*exp(stpp*(2.0*fi + bi));
        }
      }
      
	  cerr << "_"; // follow the progress with some output.

      // Fourier transform
      ffForward2.fft(ff);
      bbForward2.fft(bb);

	  cerr << ">"; // follow the progress with some output.
	
      // if this is the last slice, save the fourier transform since we just calculated it.
      if (zindex == NZ-2 && ((intstep+1) % (NZ-1) == 0))
            {
              cout << "saving farfield" << endl;
              for (ii = 0; ii < N; ii++) 
              {
                for (jj = 0; jj < N; jj++)
                {
                  x = (ii-N/2)*xstep;
                  y = (jj-N/2)*ystep;
                  if ((n = ii - N/2) < 0) n+=N;
                  if ((m = jj - N/2) < 0) m+=N;
                  farfield(ii,jj) = ff(n,m)*(1 - exp(-(x*x+y*y)/(0.55*0.55)));
                }
              }       
              sprintf(imagefile, "%03i_farfield.png", transitsCompleted);
              array_to_png(imagefile, N, farfield, (intstep+1)/(NZ-1));
              sprintf(outintfile, "%03i_ff_int.dat", transitsCompleted);
              //if (DATA) save_twod_real_array(outintfile, N, N, farfield);
              if (FFDATA) save_twod_real_array_cropped(outintfile, N, N, 103, 153, farfield);
            }

      // free-space propagate over the full step
      for (ii = 0; ii < N; ii++) 
      {
        for (jj = 0; jj < N; jj++)
        {
          ff(ii,jj) *= expDz(ii,jj);
          bb(ii,jj) *= expDz(ii,jj);
        }
      }

      // inverse Fourier transform
      ffBackward2.fftNormalized(ff);
      bbBackward2.fftNormalized(bb);
      
	  cerr << "<"; // follow the progress with some output.

      // set field values in next (previous) z-plane for fwd (bwd) field
	  #pragma omp critical
      for (ii = 0; ii < N; ii++) 
      {
        for (jj = 0; jj < N; jj++)
        {
          f_tmp(ii,jj,zindex+1) = ff(ii,jj); 
          b_tmp(ii,jj,zindex-1) = bb(ii,jj);
        }
      }

      cerr << "."; // follow the progress with some output.
    }
	// end of main parallel integration loop.

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
         sprintf(outfieldfile, "%03i_out_frame.dat", transitsCompleted);
         //sprintf(infieldfile, "%03i_in.dat", transitsCompleted);         
         if (NFDATA) save_intensity(outfieldfile, N, N, NZ-1, f);
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
  
  //array_to_png("farfield.png", N, farfield,0);
  
  /* ofstream outfile("datafile.dat");
  
  for (i = 0; i < N; i++) {
    outfile << real(f(i,N/2,NZ-1)*conj(f(i,N/2,NZ-1)))
    << "\t" << real(b(i,N/2,0)*conj(b(i,N/2,0)))
    << endl;
  } 
  */
  
  return 0;
}

#include "ImageCleaner.h"
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

#define PI	3.14159265

float *sin_x = 0;
float *cos_x = 0;
float *sin_y = 0;
float *cos_y = 0;
const int num_floats_in_cache = 4;
void cpu_fftx(float *real_image, float *imag_image, int size_x, int size_y)
{

  // Create some space for storing temporary values
  float *realOutBuffer = new float[size_y];
  float *imagOutBuffer = new float[size_y];
  // Local values
  float *fft_real = new float[size_y];
  float *fft_imag = new float[size_y];

  for(unsigned int x = 0; x < size_x; x++)
  {
    memset(realOutBuffer, 0, size_y*sizeof(float));
    memset(imagOutBuffer, 0, size_y*sizeof(float));
    for(unsigned int start_n = 0; start_n < size_y; start_n += num_floats_in_cache)
    {
      unsigned int end_n = start_n + num_floats_in_cache;
      for(unsigned int y = 0; y < size_y; y++)
      {
        for(unsigned int n = start_n; n < end_n; n++)
        {
          int term = y * n % size_y;
    realOutBuffer[y] += (real_image[x*size_y + n] * cos_y[term]) - (imag_image[x*size_y + n] * sin_y[term]);
    imagOutBuffer[y] += (imag_image[x*size_y + n] * cos_y[term]) + (real_image[x*size_y + n] * sin_y[term]);
        }
      }
    }
    // Write the buffer back to were the original values were
    memcpy(real_image+x*size_y, realOutBuffer, size_y*sizeof(float));
    memcpy(imag_image+x*size_y, imagOutBuffer, size_y*sizeof(float));
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  delete [] fft_real;
  delete [] fft_imag;
}

// This is the same as the thing above, except it has a scaling factor added to it
void cpu_ifftx(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realOutBuffer = new float[size_x];
  float *imagOutBuffer = new float[size_x];
  float *fft_real = new float[size_y];
  float *fft_imag = new float[size_y];
  for(unsigned int x = 0; x < size_x; x++)
  {
    for(unsigned int y = 0; y < size_y; y++)
    {
      for(unsigned int n = 0; n < size_y; n++)
      {
        // Compute the frequencies for this index
  float term = 2 * PI * y * n / size_y;
  fft_real[n] = cos(term);
  fft_imag[n] = sin(term);
      }

      // Compute the value for this index
      realOutBuffer[y] = 0.0f;
      imagOutBuffer[y] = 0.0f;
      for(unsigned int n = 0; n < size_y; n++)
      {
  realOutBuffer[y] += (real_image[x*size_y + n] * fft_real[n]) - (imag_image[x*size_y + n] * fft_imag[n]);
  imagOutBuffer[y] += (imag_image[x*size_y + n] * fft_real[n]) + (real_image[x*size_y + n] * fft_imag[n]);
      }

      // Incoporate the scaling factor here
      realOutBuffer[y] /= size_y;
      imagOutBuffer[y] /= size_y;
    }
    // Write the buffer back to were the original values were
    for(unsigned int y = 0; y < size_y; y++)
    {
      real_image[x*size_y + y] = realOutBuffer[y];
      imag_image[x*size_y + y] = imagOutBuffer[y];
    }
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  delete [] fft_real;
  delete [] fft_imag;
}

void cpu_ffty(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Allocate some space for temporary values
  float *realOutBuffer = new float[size_y];
  float *imagOutBuffer = new float[size_y];
  float *fft_real = new float[size_x];
  float *fft_imag = new float[size_x];
  for(unsigned int y = 0; y < size_y; y++)
  {
    for(unsigned int x = 0; x < size_x; x++)
    {
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
  float term = -2 * PI * x * n / size_x;
  fft_real[n] = cos(term);
  fft_imag[n] = sin(term);
      }

      // Compute the value for this index
      realOutBuffer[x] = 0.0f;
      imagOutBuffer[x] = 0.0f;
      for(unsigned int n = 0; n < size_x; n++)
      {
  realOutBuffer[x] += (real_image[n*size_y + y] * fft_real[n]) - (imag_image[n*size_y + y] * fft_imag[n]);
  imagOutBuffer[x] += (imag_image[n*size_y + y] * fft_real[n]) + (real_image[n*size_y + y] * fft_imag[n]);
      }
    }
    // Write the buffer back to were the original values were
    for(unsigned int x = 0; x < size_x; x++)
    {
      real_image[x*size_y + y] = realOutBuffer[x];
      imag_image[x*size_y + y] = imagOutBuffer[x];
    }
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  delete [] fft_real;
  delete [] fft_imag;
}

// This is the same as the thing about it, but it includes a scaling factor
void cpu_iffty(float *real_image, float *imag_image, int size_x, int size_y)
{
  // Create some space for storing temporary values
  float *realOutBuffer = new float[size_y];
  float *imagOutBuffer = new float[size_y];
  float *fft_real = new float[size_x];
  float *fft_imag = new float[size_x];
  for(unsigned int y = 0; y < size_y; y++)
  {
    for(unsigned int x = 0; x < size_x; x++)
    {
      // Compute the frequencies for this index
      for(unsigned int n = 0; n < size_y; n++)
      {
  // Note that the negative sign goes away for the term
  float term = 2 * PI * x * n / size_x;
  fft_real[n] = cos(term);
  fft_imag[n] = sin(term);
      }

      // Compute the value for this index
      realOutBuffer[x] = 0.0f;
      imagOutBuffer[x] = 0.0f;
      for(unsigned int n = 0; n < size_x; n++)
      {
  realOutBuffer[x] += (real_image[n*size_x + y] * fft_real[n]) - (imag_image[n*size_x + y] * fft_imag[n]);
  imagOutBuffer[x] += (imag_image[n*size_x + y] * fft_real[n]) + (real_image[n*size_x + y] * fft_imag[n]);
      }

      // Incorporate the scaling factor here
      realOutBuffer[x] /= size_x;
      imagOutBuffer[x] /= size_x;
    }
    // Write the buffer back to were the original values were
    for(unsigned int x = 0; x < size_x; x++)
    {
      real_image[x*size_y + y] = realOutBuffer[x];
      imag_image[x*size_y + y] = imagOutBuffer[x];
    }
  }
  // Reclaim some memory
  delete [] realOutBuffer;
  delete [] imagOutBuffer;
  delete [] fft_real;
  delete [] fft_imag;
}

void cpu_filter(float *real_image, float *imag_image, int size_x, int size_y)
{
  int eightX = size_x/8;
  int eight7X = size_x - eightX;
  int eightY = size_y/8;
  int eight7Y = size_y - eightY;
  for(unsigned int x = 0; x < size_x; x++)
  {
    for(unsigned int y = 0; y < size_y; y++)
    {
      if(!(x < eightX && y < eightY) &&
   !(x < eightX && y >= eight7Y) &&
   !(x >= eight7X && y < eightY) &&
   !(x >= eight7X && y >= eight7Y))
      {
  // Zero out these values
  real_image[y*size_x + x] = 0;
  imag_image[y*size_x + x] = 0;
      }
    }
  }
}

float transpose(float *image, float *image_transpose, int size_x, int size_y)
{
  for(int x = 0; x < size_x; x++){
    for (int y = 0; y < size_y; y++){
      image_transpose[y*size_x+x] = image[x*size_y+y];
    }
  }
}
float mimageCleaner(float *real_image, float *imag_image, int size_x, int size_y)
{
  
  sin_x = new float[size_x];
  cos_x = new float[size_x];
  unsigned int size = size_x * size_y;
  float *real_image_transpose = new float[size];
  float *imag_image_transpose = new float[size];
  int eightX = size_x/8;
  int eight7X = size_x - eightX;
  int eightY = size_y/8;
  int eight7Y = size_y - eightY;
  #pragma omp parallel
  {
    // calculate sin_x, sin_y, cos_x, cos_y
    #pragma omp for
    for(int i = 0; i < size_x; i++){
      float term = -2.0 * PI * i / size_x;
      sin_x[i] = sin(term);
      cos_x[i] = cos(term);
    }
    if (size_y == size_x){
      sin_y = sin_x;
      cos_y = cos_x;
    }
    else{
      sin_y = new float[size_y];
      cos_y = new float[size_y];
      #pragma omp for
      for (int i = 0; i < size_y; i++){
        float term = -2.0 * PI * i / size_y;
        sin_y[i] = sin(term);
        cos_y[i] = cos(term);
      }
    }

    // cpu_fftx
    #pragma omp for
    for(unsigned int x = 0; x < size_x; x++)
    {
      for(unsigned int y = 0; y < size_y; y++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        // #pragma omp for         
        // #pragma omp for nowait reduction(+:real) reduction(+:imag) 
        for(unsigned int n = 0; n < size_y; n++)
        {
          int term = y * n % size_y;
          int idx = x*size_y + n;
          real = real + (real_image[idx] * cos_y[term]) - (imag_image[idx] * sin_y[term]);
          imag = imag + (imag_image[idx] * cos_y[term]) + (real_image[idx] * sin_y[term]);
        }
        real_image_transpose[y*size_x + x] = real;
        imag_image_transpose[y*size_x + x] = imag;
      }
    }

    // cpu_ffty
    #pragma omp for nowait
    for(unsigned int y = 0; y < eightY; y++)
    {
      for(unsigned int x = 0; x < eightX; x++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        for(unsigned int n = 0; n < size_x; n++)
        {
          int term = x * n % size_x;
          int idx = y * size_x + n;
          real += (real_image_transpose[idx] * cos_x[term]) - (imag_image_transpose[idx] * sin_x[term]);
          imag += (imag_image_transpose[idx] * cos_x[term]) + (real_image_transpose[idx] * sin_x[term]);
        }
        real_image[x*size_y+y] = real;
        imag_image[x*size_y+y] = imag;
      }
      for(unsigned int x = eight7X; x < size_x; x++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        for(unsigned int n = 0; n < size_x; n++)
        {
          int term = x * n % size_x;
          int idx = y * size_x + n;
          real += (real_image_transpose[idx] * cos_x[term]) - (imag_image_transpose[idx] * sin_x[term]);
          imag += (imag_image_transpose[idx] * cos_x[term]) + (real_image_transpose[idx] * sin_x[term]);
        }
        real_image[x*size_y+y] = real;
        imag_image[x*size_y+y] = imag;
      }
    }

    #pragma omp for
    for(unsigned int y = eight7Y; y < size_y; y++)
    {
      for(unsigned int x = 0; x < eightX; x++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        for(unsigned int n = 0; n < size_x; n++)
        {
          int term = x * n % size_x;
          int idx = y * size_x + n;
          real += (real_image_transpose[idx] * cos_x[term]) - (imag_image_transpose[idx] * sin_x[term]);
          imag += (imag_image_transpose[idx] * cos_x[term]) + (real_image_transpose[idx] * sin_x[term]);
        }
        real_image[x*size_y+y] = real;
        imag_image[x*size_y+y] = imag;
      }
      for(unsigned int x = eight7X; x < size_x; x++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        for(unsigned int n = 0; n < size_x; n++)
        {
          int term = x * n % size_x;
          int idx = y * size_x + n;
          real += (real_image_transpose[idx] * cos_x[term]) - (imag_image_transpose[idx] * sin_x[term]);
          imag += (imag_image_transpose[idx] * cos_x[term]) + (real_image_transpose[idx] * sin_x[term]);
        }
        real_image[x*size_y+y] = real;
        imag_image[x*size_y+y] = imag;
      }
    }

    // cpu_iffx
    #pragma omp for nowait
    for(unsigned int x = 0; x < eightX; x++)
    {
      for(unsigned int y = 0; y < size_y; y++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        // #pragma omp for         
        // #pragma omp for nowait reduction(+:real) reduction(+:imag) 
        for(unsigned int n = 0; n < eightY; n++)
        {
          int term = y * n % size_y;
          int idx = x*size_y + n;
          real = real + (real_image[idx] * cos_y[term]) + (imag_image[idx] * sin_y[term]);
          imag = imag + (imag_image[idx] * cos_y[term]) - (real_image[idx] * sin_y[term]);
        }
        for(unsigned int n = eight7Y; n < size_y; n++)
        {
          int term = y * n % size_y;
          int idx = x*size_y + n;
          real = real + (real_image[idx] * cos_y[term]) + (imag_image[idx] * sin_y[term]);
          imag = imag + (imag_image[idx] * cos_y[term]) - (real_image[idx] * sin_y[term]);
        }
        real_image_transpose[y*size_x + x] = real/size_y;
        imag_image_transpose[y*size_x + x] = imag/size_y;
      }
    }
    #pragma omp for
    for(unsigned int x = eight7X; x < size_x; x++)
    {
      for(unsigned int y = 0; y < size_y; y++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        // #pragma omp for         
        // #pragma omp for nowait reduction(+:real) reduction(+:imag) 
        for(unsigned int n = 0; n < eightY; n++)
        {
          int term = y * n % size_y;
          int idx = x*size_y + n;
          real = real + (real_image[idx] * cos_y[term]) + (imag_image[idx] * sin_y[term]);
          imag = imag + (imag_image[idx] * cos_y[term]) - (real_image[idx] * sin_y[term]);
        }
        for(unsigned int n = eight7Y; n < size_y; n++)
        {
          int term = y * n % size_y;
          int idx = x*size_y + n;
          real = real + (real_image[idx] * cos_y[term]) + (imag_image[idx] * sin_y[term]);
          imag = imag + (imag_image[idx] * cos_y[term]) - (real_image[idx] * sin_y[term]);
        }
        real_image_transpose[y*size_x + x] = real/size_y;
        imag_image_transpose[y*size_x + x] = imag/size_y;
      }
    }

    // cpu_iffy
    #pragma omp for
    for(unsigned int y = 0; y < size_y; y++)
    {
      for(unsigned int x = 0; x < size_x; x++)
      {
        float real = 0.0f;
        float imag = 0.0f;
        for(unsigned int n = 0; n < eightX; n++)
        {
          int term = x * n % size_x;
          int idx = y * size_x + n;
          real += (real_image_transpose[idx] * cos_x[term]) + (imag_image_transpose[idx] * sin_x[term]);
          imag += (imag_image_transpose[idx] * cos_x[term]) - (real_image_transpose[idx] * sin_x[term]);
        }
        for(unsigned int n = eight7X; n < size_x; n++)
        {
          int term = x * n % size_x;
          int idx = y * size_x + n;
          real += (real_image_transpose[idx] * cos_x[term]) + (imag_image_transpose[idx] * sin_x[term]);
          imag += (imag_image_transpose[idx] * cos_x[term]) - (real_image_transpose[idx] * sin_x[term]);
        }
        real_image[x*size_y+y] = real/size_x;
        imag_image[x*size_y+y] = imag/size_x;
      }
    }
  }
}
float imageCleaner(float *real_image, float *imag_image, int size_x, int size_y)
{
  // These are used for timing
  struct timeval tv1, tv2;
  struct timezone tz1, tz2;


  sin_x = new float[size_x];
  cos_x = new float[size_x];
  for(int i = 0; i < size_x; i++){
    float term = -2.0 * PI * i / size_x;
    sin_x[i] = sin(term);
    cos_x[i] = cos(term);
  }
  if (size_y == size_x){
    sin_y = sin_x;
    cos_y = cos_x;
  }
  else{
    sin_y = new float[size_y];
    cos_y = new float[size_y];
    #pragma omp for
    for (int i = 0; i < size_y; i++){
      float term = -2.0 * PI * i / size_y;
      sin_y[i] = sin(term);
      cos_y[i] = cos(term);
    }
  }
  // Start timing
  gettimeofday(&tv1,&tz1);

  // mimageCleaner(real_image, imag_image, size_x, size_y);
  // Perform fft with respect to the x direction
  cpu_fftx(real_image, imag_image, size_x, size_y);
  // Perform fft with respect to the y direction
  // cpu_ffty(real_image, imag_image, size_x, size_y);

  // Filter the transformed image
  // cpu_filter(real_image, imag_image, size_x, size_y);

  // Perform an inverse fft with respect to the x direction
  // cpu_ifftx(real_image, imag_image, size_x, size_y);
  // Perform an inverse fft with respect to the y direction
  // cpu_iffty(real_image, imag_image, size_x, size_y);

  // End timing
  gettimeofday(&tv2,&tz2);

  // Compute the time difference in micro-seconds
  float execution = ((tv2.tv_sec-tv1.tv_sec)*1000000+(tv2.tv_usec-tv1.tv_usec));
  // Convert to milli-seconds
  execution /= 1000;
  // Print some output
  printf("OPTIMIZED IMPLEMENTATION STATISTICS:\n");
  printf("  Optimized Kernel Execution Time: %f ms\n\n", execution);
  return execution;
}

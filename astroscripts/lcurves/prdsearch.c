#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdbool.h>

typedef enum {
    LSPOWER
} KnownRoutines;

typedef struct s_thrargs {
    KnownRoutines rname;   // Routine name
    float *inX;
    float *inY;
    float *outX;
    float *outY;
    uint inlen;            // Length of the input arrays
    uint outlen;           // Length of the output arrays 
    bool progress;
    uint res;              // Store the value returned by routine
} thrargs;   // Arguments to pass a thread routine



/*Calculate raw (non-normalized) Lomb-Scargle power at
arbitrary frequencies. See Horne & Baliunas, 1986, 
%ApJ, 302, 757
http://adsabs.harvard.edu/abs/1986ApJ...302..757H
The arguments are:
  time  - time stamps of the light curve;
  data  - count rate or flux;
  freq  - desired frequencies;
  power - an array to write the result, must have
          the same length as freq;
  inlen - length of the time and data arrays; 
  oulen - length of the freq and power arrays;
  progress - it true, prints the progress.

Returns the number of processed frequences.
*/
static uint _lspower(float *time, float *data, float *freq, 
            float *spec, uint inlen, uint oulen, bool progress)
{
    uint i, j;
    float sins, coss;  // sum of sin and cos in eq. (2) 
    float omtau;       // omega*tau in eq. (2)
    float arg;         // omega*(t-tau)
    float sina, ysin, sinsq, cosa, ycos, cossq;

    // printf("I'm in lspower %d %f\n", oulen, freq[0]);
    // fflush(stdout);

    //Number of frequencies to treat to change progress by 0.01%
    uint onemark = ceilf((float)oulen/10000);
    if (progress == true)  // Print progress
    {
        printf("Progress: 00.00%%");
        fflush(stdout);
    }

    for (i=0; i<oulen; i++)   // Frequencies
    {
        sins = coss = omtau = 0.0f;
        for (j=0; j<inlen; j++)   // Calculate phase shift
        {
            sins += sinf(4.0f*M_PI*freq[i]*time[j]);
            coss += cosf(4.0f*M_PI*freq[i]*time[j]);
        }
        omtau = 0.5f*atan2f(sins, coss);

        ycos = ysin = cossq = sinsq = 0.0f;
        for (j=0; j<inlen; j++)
        {
            arg = 2.0f*M_PI*freq[i]*time[j] - omtau;
            cosa = cosf(arg);
            sina = sinf(arg);
            ycos += data[j]*cosa;
            ysin += data[j]*sina;
            cossq += powf(cosa, 2);
            sinsq += powf(sina, 2);
        }
        spec[i] = powf(ycos, 2)/cossq +
            powf(ysin, 2)/sinsq;

        if (progress == true)
            if (i % onemark == 0)
            {
                printf("\b\b\b\b\b\b%5.2f%%", (100.0f*i)/oulen);
                fflush(stdout);
            }
    }
    if (progress == true)
    {
        printf("\b\b\b\b\b\b%5.2f%%\n", 100.0f);
        fflush(stdout);
    }
    return i;
}


static void* _wrapper(void *arguments)
{
    thrargs *args = (thrargs*) arguments;
    switch(args->rname)
    {
        case LSPOWER: 
            args->res = _lspower(args->inX, args->inY, args->outX,
                args->outY, args->inlen, args->outlen, args->progress);
    }
    pthread_exit(NULL); 
}


/*Provide multi-threading for routines.*/
static uint _threader(KnownRoutines rname, float *inX, float *inY,
    float *outX, float *outY, uint inlen, uint outlen, bool progress, 
    unsigned char nth)
{

    uint i, fullres = 0, chunksize = outlen / nth;
    pthread_t thid[nth];     // Thread ID
    thrargs args[nth];    // Arguments for each worker

    for (i=0; i<nth; i++)  // Create threads
    {
        //Prepare arguments
        args[i].rname = rname;
        args[i].inX = inX;
        args[i].inY = inY;
        args[i].inlen = inlen;
        args[i].outX = &outX[i*chunksize];
        args[i].outY = &outY[i*chunksize];
        args[i].outlen = chunksize;
        args[i].progress = false;
        if (i == (nth-1))  // this is the ultimate runner
        {
            args[i].outlen = outlen - i*chunksize;  // Get all the remainig points
            args[i].progress = progress;            // and turn on the progress
        }
            
        if (pthread_create(&thid[i], NULL, _wrapper, &args[i]))
            return 0;
    }

    for (i=0; i<nth; i++)  // Collect the results
    {
        pthread_join(thid[i], NULL);
        fullres += args[i].res;
    }
    return fullres;
}


uint lspower(float *time, float *data, float *freq, 
            float *spec, uint inlen, uint oulen,
            bool progress, unsigned char nth)
{

    if (nth == 1)   // We don't need threading
        return _lspower(time, data, freq, spec, inlen,
             oulen, progress);
    else
        return _threader(LSPOWER, time, data, freq, spec, inlen,
             oulen, progress, nth);

}
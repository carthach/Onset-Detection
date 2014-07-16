#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <fftw3.h>
#include <sndfile.h>

#define MAX( a, b ) ( ( a > b) ? a : b )
#define MIN( a, b ) ( ( a < b) ? a : b )

#define	BLOCK_SIZE 256

//PARAM
float threshold;

//FFTW
fftw_complex	*data, *fft_result, *ifft_result;
fftw_plan	plan_forward, plan_backward;
int		i;

//libsndfile
SNDFILE  *infile = NULL ;
SNDFILE	 *outwav = NULL ;
SF_INFO	sfinfo ;

char *progname, *infilename, *outfilename ;
FILE *outfile = NULL ;


FILE *onsetTimesFile;


float last_fft[BLOCK_SIZE/2];
int noOfFrames;
int averageWindowSize;
int sizeOfFunctions;

float *onsetFunction, *averageFunction;

static void print_usage (char *progname)
{	
    printf ("\nUsage : %s <input file> <output file>\n", progname) ;
    puts ("\n"
	  "    Where the output file will contain a line for each frame\n"
	  "    and a column for each channel.\n"
	) ;

} /* print_usage */

/*convert to mono*/
float f_mono(float l, float r){
    return(l+r)/2.0f;
}

/*Euclidean distance between frames*/
float distance(float *array0, float *array1, int array_size)
{
    float dist;
    int i;
    for(i=0; i<array_size; i++)
    {
	float val = array0[i]-array1[i];
	val = val*val;
	dist += val;
    }
    dist = sqrt(dist);

    return dist;
}

/*Euclidean distance between frames*/
float difference(float *array0, float *array1, int array_size)
{
    float dist;
    int i;
    for(i=0; i<array_size; i++)
    {
	float val = array0[i]-array1[i];
	dist += val;
    }

    return dist;
}

float rollingAverage()
{
    int i;

    for(i=0; i<sizeOfFunctions;i++)
    {
	int start = MAX(0, i-averageWindowSize);
	int end = MIN(sizeOfFunctions-1, i+averageWindowSize);

	int j;
	float mean;
	for(j=start; j<end; j++){
	    mean += onsetFunction[j];
	}
	mean /= (float)(end-start);

	averageFunction[i] = mean * threshold;

    }
}

float removeThreshold()
{
    int i;
    for(i=0; i<sizeOfFunctions;i++) {
	if(onsetFunction[i] < averageFunction[i]) {
	    onsetFunction[i] = 0.0;
	}
    }	
}

float peakDetection()
{
    int i,j;

    int w=20;
    
    for(i=0; i<sizeOfFunctions-1; i++) {
	int start = MAX(0, i-w);
	int end = MIN(sizeOfFunctions-1, i+w);

	if(onsetFunction[i] > averageFunction[i]) {
	    for(j=start; j<end; j++) {
		if(onsetFunction[i]<onsetFunction[j]) {
		    onsetFunction[i] = 0.0f;
		    break;
		}
	    }
	} else {
	    onsetFunction[i] = 0.0f;
	}
    }
}

void printArray(FILE* fileToPrint, float* arrayToPrint, int arraySize)
{
    int i=0;
    for(i=0; i<arraySize;i++)
	fprintf(fileToPrint, "%f\n", arrayToPrint[i]);
}


static void process ()
{
    float buf [sfinfo.channels * BLOCK_SIZE];
    float outbuf[BLOCK_SIZE];
    float current_fft[BLOCK_SIZE/2];
    int k, m, readcount ;

    //Get every frame
    while ((readcount = sf_readf_float (infile, buf, BLOCK_SIZE)) > 0)
    {
	for (k = 0 ; k < readcount ; k++)
	{	for (m = 0 ; m < sfinfo.channels ; m++){ 

		float mono_sample = buf[k];
		if(sfinfo.channels==2)
		    mono_sample = f_mono(buf[k*2],buf[k*2+1]);

		float window_multiplier = 0.5 * (1.0-cos(2.0*M_PI*(float)k/(float)BLOCK_SIZE));

		outbuf[k] = mono_sample * window_multiplier;

		data[k][0] = outbuf[k];
		data[k][1] = 0.0;
	    }
	}

	fftw_execute(plan_forward);
	int f;

	for(f=0; f<BLOCK_SIZE/2; f++){
	    float mag = (float)sqrt(fft_result[f][0]*fft_result[f][0]+fft_result[f][1]*fft_result[f][1]);
	    current_fft[f] = mag;
	}

//	float dist = distance(current_fft, last_fft, SIZE/2 );
	float dist = difference(current_fft, last_fft, BLOCK_SIZE/2 );

	onsetFunction[noOfFrames] = dist;
	
	memcpy(current_fft, last_fft, sizeof(current_fft));
	sf_writef_float(outwav, outbuf, BLOCK_SIZE);
	noOfFrames++;
    }

    rollingAverage();
//    removeThreshold();
    peakDetection();
//    printArray(onsetFunctionFile, onsetFunction, sizeOfFunctions);

    int i;
    for(i=0; i<sizeOfFunctions; i++)
	if(onsetFunction[i] > 0.0)
	    fprintf(outfile, "%f\n", i * (float)BLOCK_SIZE/(float)sfinfo.samplerate);

    return ;
} /* convert_to_text */

int main (int argc, char * argv [])
{



    progname = strrchr (argv [0], '/') ;
    progname = progname ? progname + 1 : argv [0] ;

    if (argc < 3)
    {	print_usage (progname) ;
	return 1 ;
    } ;

    infilename = argv [1] ;
    outfilename = argv [2] ;

    if (strcmp (infilename, outfilename) == 0)
    {	printf ("Error : Input and output filenames are the same.\n\n") ;
	print_usage (progname) ;
	return 1 ;
    } ;

    if (infilename [0] == '-')
    {	printf ("Error : Input filename (%s) looks like an option.\n\n", infilename) ;
	print_usage (progname) ;
	return 1 ;
    } ;

    if (outfilename [0] == '-')
    {	printf ("Error : Output filename (%s) looks like an option.\n\n", outfilename) ;
	print_usage (progname) ;
	return 1 ;
    } ;

    if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
    {	printf ("Not able to open input file %s.\n", infilename) ;
	puts (sf_strerror (NULL)) ;
	return 1 ;
    } ;

    /* Open the output file. */
    if ((outfile = fopen (outfilename, "w")) == NULL)
    {	printf ("Not able to open output file %s : %s\n", outfilename, sf_strerror (NULL)) ;
	return 1 ;
    } ;

    if (argc == 4)
	threshold = atof(argv[3]);
    else
	threshold = 1.5;

    noOfFrames = 0;

    data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * BLOCK_SIZE);
    fft_result  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * BLOCK_SIZE);
    ifft_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * BLOCK_SIZE);
 
    plan_forward  = fftw_plan_dft_1d(BLOCK_SIZE, data, fft_result,
				     FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(BLOCK_SIZE, fft_result, ifft_result,
				     FFTW_BACKWARD, FFTW_ESTIMATE);

    SF_INFO outsfinfo = sfinfo;
    outsfinfo.channels = 1;



    float calc = (float)sfinfo.frames/(float)BLOCK_SIZE;
    calc = ceil(calc);
    printf("calc is %f\n", calc);
    sizeOfFunctions = calc;
    printf("sizeOfFunctions is %d\n", sizeOfFunctions);


    onsetFunction = malloc(sizeof(float)*sizeOfFunctions);
    averageFunction = malloc(sizeof(float)*sizeOfFunctions);

    averageWindowSize = 43;


    process() ;

    printf("frames: %d\n", sfinfo.frames);
    printf("blocks: %f\n", ceil((float)sfinfo.frames/(float)BLOCK_SIZE));
    printf("noOfFrames: %d\n", noOfFrames);
    printf("averageWindowSize: %d\n", averageWindowSize);

    sf_close(infile);
    fclose(outfile);

    FILE *onsetFunctionFile;
    FILE *averageFunctionFile;

    onsetFunctionFile = fopen("onsetFunction.txt", "w");
    averageFunctionFile = fopen("averageFunction.txt", "w");
    printArray(onsetFunctionFile, onsetFunction, sizeOfFunctions);
    printArray(averageFunctionFile, averageFunction, sizeOfFunctions);
    fclose (onsetFunctionFile);
    fclose (averageFunctionFile);

    free(onsetFunction);
    free(averageFunction);
    return 0 ;
} /* main */


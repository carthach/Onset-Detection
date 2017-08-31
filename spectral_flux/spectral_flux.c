#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include <sndfile.h>

#define MAX( a, b ) ( ( a > b) ? a : b )
#define MIN( a, b ) ( ( a < b) ? a : b )

#define DEFAULT_THRESHOLD 1.5
#define	DEFAULT_BLOCK_SIZE 256
#define	DEFAULT_AVERAGE_WINDOW_SIZE 43

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


FILE *onset_times_file;

int block_size;
int block_count;
int average_window_size;
int function_size;

float *onset_function, *average_function, *last_fft;

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
    float dist = 0.0;
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
    float dist = 0.0;
    int i;
    for(i=0; i<array_size; i++)
    {
        float val = array0[i]-array1[i];
        dist += val;
    }
    return dist;
}

void rollingAverage()
{
    int i;
    
    for(i=0; i<function_size;i++)
    {
        int start = MAX(0, i-average_window_size);
        int end = MIN(function_size-1, i+average_window_size);
        
        int j;
        float mean = 0.0;
        for(j=start; j<end; j++){
            mean += onset_function[j];
        }
        mean /= (float)(end-start);
        
        average_function[i] = mean * threshold;
    }
}

void removeThreshold()
{
    int i;
    for(i=0; i<function_size;i++) {
        if(onset_function[i] < average_function[i]) {
            onset_function[i] = 0.0;
        }
    }
}

void peakDetection()
{
    int i,j;
    
    int w=20;
    
    for(i=0; i<function_size-1; i++) {
        int start = MAX(0, i-w);
        int end = MIN(function_size-1, i+w);
        
        if(onset_function[i] > average_function[i]) {
            for(j=start; j<end; j++) {
                if(onset_function[i]<onset_function[j]) {
                    onset_function[i] = 0.0f;
                    break;
                }
            }
        } else {
            onset_function[i] = 0.0f;
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
    float buf [sfinfo.channels * block_size];
    float outbuf[block_size];
    float current_fft[block_size/2];
    int k, m;
    long readcount;
    
    //Get every frame
    while ((readcount = sf_readf_float (infile, buf, block_size)) > 0)
    {
        for (k = 0 ; k < readcount ; k++)
        {
            for (m = 0 ; m < sfinfo.channels ; m++)
            {
                //Convert to mono
                float mono_sample = buf[k];
                
                if(sfinfo.channels==2)
                    mono_sample = f_mono(buf[k*2],buf[k*2+1]);
                
                //Window
                float window_multiplier = 0.5 * (1.0-cos(2.0*M_PI*(float)k/(float)block_size));
    
                outbuf[k] = mono_sample * window_multiplier;
                
                //Eh....
                data[k][0] = outbuf[k];
                data[k][1] = 0.0;
            }
        }
        
        //Get FFT
        fftw_execute(plan_forward);
        int f;
        
        //Retain Magnitudes
        for(f=0; f<block_size/2; f++){
            float mag = (float)sqrt(fft_result[f][0]*fft_result[f][0]+fft_result[f][1]*fft_result[f][1]);
            current_fft[f] = mag;
        }
        
        //	float dist = distance(current_fft, last_fft, SIZE/2 );
        float dist = difference(current_fft, last_fft, block_size/2);
        
        onset_function[block_count] = dist;
        
        memcpy(current_fft, last_fft, sizeof(current_fft));
        sf_writef_float(outwav, outbuf, block_size);
        block_count++;
    }
    
    rollingAverage();
    //    removeThreshold();
    peakDetection();
    //    printArray(onset_functionFile, onset_function, function_size);
    
    int i;
    for(i=0; i<function_size; i++)
        if(onset_function[i] > 0.0)
            fprintf(outfile, "%f\n", i * (float)block_size/(float)sfinfo.samplerate);
    
    return ;
} /* convert_to_text */

int main (int argc, char * argv [])
{
    //Set defaults
    threshold = DEFAULT_THRESHOLD;
    block_size = DEFAULT_BLOCK_SIZE;
    average_window_size = DEFAULT_AVERAGE_WINDOW_SIZE;
    block_count = 0;
    
    //Don't know what this is ha
    progname = strrchr (argv [0], '/') ;
    progname = progname ? progname + 1 : argv [0] ;
    
    //Explain instructions
    if (argc < 3)
    {
        print_usage (progname) ;
        return 1 ;
    };
    
    /* Process Optional Arguments */
    
    int index;
    int c;
    int errflg = 0;
    
    while ((c = getopt (argc, argv, "b:t:")) != -1)
        switch (c)
    {
        case 'b':
            block_size = atoi(optarg);
            printf("Block Size: %d\n", block_size);
            break;
        case 'w':
            average_window_size = atoi(optarg);
            printf("Moving Average Window Size: %d\n", average_window_size);
            break;
        case 't':
            threshold = atof(optarg);
            printf("Threshold: %f\n", threshold);
            break;
        case ':':       /* -f or -o without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr,
                    "Unrecognized option: '-%c'\n", optopt);
            errflg++;
            return 1;
        default:
            abort ();
    }
    
    /* Process the file arguments */
    
    for (index = optind; index < argc; index++)
        printf ("Non-option argument %s\n", argv[index]);
    
    //Get the filenames
    infilename = argv [optind] ;
    outfilename = argv [optind+1] ;
    
    //Some unnecessary file checking
    if (strcmp (infilename, outfilename) == 0)
    {
        printf ("Error : Input and output filenames are the same.\n\n") ;
        print_usage (progname) ;
        return 1 ;
    };
    
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
    
    /* Get file handles */
    
    //Open the audio file
    if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
    {	printf ("Not able to open input file %s.\n", infilename) ;
        puts (sf_strerror (NULL)) ;
        return 1 ;
    } ;
    
    //Open the outputfile
    if ((outfile = fopen (outfilename, "w")) == NULL)
    {	printf ("Not able to open output file %s : %s\n", outfilename, sf_strerror (NULL)) ;
        return 1 ;
    } ;
    
    /* FFT Allocation */
    
    last_fft = malloc(sizeof(float)*block_size);
    
    data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * block_size);
    fft_result  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * block_size);
    ifft_result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * block_size);
    
    plan_forward  = fftw_plan_dft_1d(block_size, data, fft_result,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(block_size, fft_result, ifft_result,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    
    SF_INFO outsfinfo = sfinfo;
    outsfinfo.channels = 1;
    
    //Get the number of blocks we need process for the file
    float blocks_per_file = (float)sfinfo.frames/(float)block_size;
    printf("blocks_per_file: %f\n", blocks_per_file);
    
    //This then becomes the function size
    function_size = (int)ceil(blocks_per_file);
    printf("function_size: %d\n", function_size);
    
    //Allocate
    onset_function = malloc(sizeof(float)*function_size);
    average_function = malloc(sizeof(float)*function_size);
    
    //Get to work!
    process() ;
    
    printf("frames: %d\n", (int)sfinfo.frames);
    printf("blocks_processed: %d\n", block_count);
    printf("average_window_size: %d\n", average_window_size);
    
    /* Cleanup */
    sf_close(infile);
    fclose(outfile);
    
    //Debugging the various functions
    FILE *onset_functionFile = fopen("onset_function.txt", "w");
    FILE *average_functionFile = fopen("average_function.txt", "w");
    printArray(onset_functionFile, onset_function, function_size);
    printArray(average_functionFile, average_function, function_size);
    fclose (onset_functionFile);
    fclose (average_functionFile);
    
    free(onset_function);
    free(average_function);
    free(last_fft);
    return 0 ;
} /* main */


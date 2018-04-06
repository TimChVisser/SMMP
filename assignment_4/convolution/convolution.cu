#include <stdio.h>
#include <string.h>
#include "timer.h"

#define image_height 8192
#define image_width 1024
#define filter_height 5
#define filter_width 5
#define filter_width_sep 3

#define border 2
#define border_height 4
#define border_width 4
#define input_height (image_height + border_height)
#define input_width (image_width + border_width)
#define input_height_sep (image_height + 2)
#define input_width_sep (image_width + 2)

#define block_size_x 32
#define block_size_y 32

#define block_size_hor_x 32
#define block_size_hor_y 8

#define block_size_ver_x 8
#define block_size_hor_y 32
// #define DEBUG
#define SEED 1234


#define block_size_x1 16
#define block_size_y1 16

using namespace std;

__constant__ float d_kernel[filter_width*filter_height];

void convolutionSeq(float *output, float *input, float *filter) {
    //for each pixel in the output image

  timer sequentialTime = timer("Sequential");

  sequentialTime.start();

    for (int y=0; y < image_height; y++) {
        for (int x=0; x < image_width; x++) {

            //for each filter weight
            for (int i=0; i < filter_height; i++) {
                for (int j=0; j < filter_width; j++) {
                    output[y*image_width+x] += input[(y+i)*input_width+x+j]
                        * filter[i*filter_width+j];
                }
            }
            output[y*image_width+x] /= 35;
        }
    }
  sequentialTime.stop();
  cout << "convolution (sequential): \t\t" << sequentialTime << endl;

}


__global__ void convolution_kernel_naive(float *output, float *input, float * filter) {
    // global mem address for this thread
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int x = blockIdx.x * blockDim.x + threadIdx.x;

    // if(x>input_width-1 || y > input_height -1){
    //     return;
    // }
    //fprintf("error %i:%i",y,x);
    float sum = 0;
    //for each filter weight
    for (int i=0; i < filter_height; i++) {
        for (int j=0; j < filter_width; j++) {
           sum += input[(y+i)*input_width+x+j] * filter[i*filter_width+j];
        }
    }
    output[y*image_width+x] = sum/35;

}

__global__ void convolution_kernel_shared(float *output, float *input) {
    // global mem address for this thread
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    // if(x<image_width-1 && y < image_height-1){
        const size_t shared_width = block_size_x1+border_width;
        const size_t shared_height= block_size_y1+ border_height;

        __shared__ float data[shared_width * shared_height];
        //loading top left tile
        data[threadIdx.y*shared_width+threadIdx.x] = input[y*input_width+x];
        // // loading boottom left tile
        if(threadIdx.y < border_width){
            data[(threadIdx.y+blockDim.y)*shared_width+threadIdx.x]
                                    = input[(y+blockDim.y)*input_width+x];
        }
        // loading upper right
        if(threadIdx.x < border_height){
            data[(threadIdx.y)*shared_width+threadIdx.x+blockDim.x]
                                        = input[y*input_width+x+blockDim.x];
        }
        // loading bottom right
        if(threadIdx.x < border_width && threadIdx.y < border_height){
            data[(threadIdx.y+blockDim.y)*shared_width+threadIdx.x+blockDim.x]
                                        = input[(y+blockDim.y)*input_width+x+blockDim.x];
        }

        __syncthreads();
        //fprintf("error %i:%i",y,x);
        float sum = 0;
        //for each filter weight
        for (int i=0; i < filter_height; i++) {
            for (int j=0; j < filter_width; j++) {
            sum += data[(threadIdx.y+i)*shared_width+threadIdx.x+j] * d_kernel[i*filter_width+j];
            }
        }
        output[y*image_width+x] = sum/35;
    // }
}




// __global__ void convolution_kernel_horizontal(float *output, float *input) {
//     // global mem address for this thread
//     const int y = blockIdx.y * blockDim.y + threadIdx.y;
//     const int x = blockIdx.x * blockDim.x + threadIdx.x;
//     const size_t kernel_half = 2;
//     const size_t kernel_size = 5;
//     const size_t width = block_size_hor_x + kernel_size;
//     __shared__ float data[block_size_hor_y * width];
//      //loading top left row tile
//     data[threadIdx.x] = input[y*input_width+x];
//     //loading top right row tile
//     if(threadIdx.x < kernel_size){
//         data[threadIdx.x+blockDim.x] = input[y*input_width+x + blockDim.x];
//     }
//     __syncthreads();
//     float sum = 0;
//     //for each filter weight
//     int start = y*input_width_sep+x;
//     for (int j=0; j < kernel_size; j++) {
//         sum += data[j];

//     }
//     output[y*input_width_sep+x] = sum;

// }

// __global__ void convolution_kernel_vertical(float *output, float *input) {
//     // global mem address for this thread
//     const int y = blockIdx.y * blockDim.y + threadIdx.y;
//     const int x = blockIdx.x * blockDim.x + threadIdx.x;
//     const size_t kernel_half = 2;
//     const size_t kernel_size = 5;

//     if(x>input_width_sep-1 || y > input_height -1){
//         return;
//     }
//      __shared__ float data[blockDim.y + kernel_size];
//     //loading top column tile
//     data[threadIdx.y] = input[y*input_width+x];
//     //loading botom column tile
//     if(threadIdx.x < kernel_size){
//         data[threadIdx.y+blockDim.y] = input[y*input_width+x + blockDim.x];
//     }
//     float sum = 0;
//     //for each filter weight
//     for (int i=0; i < filter_width_sep; i++) {
//         sum += input[(y+i)*input_width_sep+x] * d_kernelS[i];
//     }
//     output[y*image_width+x] = sum;

// }

// __global__ void convolution_kernel_horizontal_shared(float *output, float *input) {
//     // global mem address for this thread
//     const int y = blockIdx.y * blockDim.y + threadIdx.y;
//     const int x = blockIdx.x * blockDim.x + threadIdx.x;

//     __shared__ float data[block_size_hor+2];
//     data[blockIdx.x] = input[y*input_width+x];
//     if(blockIdx.x > blockDim.x -2)
//     {
//         int stride = blockDim.x - blockIdx.x
//         data[blockIdx.x +  stride] =  input[y*input_width+x + stride];
//     }
//     __syncthreads();

//     float sum = 0;
//     //for each filter weight
//     if(blockIdx.x < blockDim.x -1 || x > inpu_width)
//     int i = 0;
//     for (int j=0; j < filter_width; j++) {
//         sum += data[(y+i)*input_width+x+j] * d_kernelS[j];
//     }
//     output[y*input_width+x] = sum;

// }

// __device__ void
// d_boxfilter_x(float * input, float *output, int width, int height, int radius)
// {
//     float scale = 1.0f / (float)((r * 2) + 1);

//     float t;
//     // do left edge
//     t = input[0] * r;

//     for (int x = 0; x < w ; x++)
//     {
//         t += id[x + r];
//         t -= id[x - r - 1];
//         od[x] = t * scale;
//     }

// }

// // process column
// __device__ void
// d_boxfilter_y(float *id, float *od, int w, int h, int r)
// {
//     const int y = blockIdx.y * blockDim.y + threadIdx.y;
//     const int x = blockIdx.x * blockDim.x + threadIdx.x;

//     float scale = 1.0f / (float)((r << 1) + 1);

//     float t;
//     // do left edge
//     t = id[0] * r;

//     for (int y = 0; y < (r + 1); y++)
//     {
//         t += id[y * w];
//     }

//     od[0] = t * scale;

//     for (int y = 1; y < (r + 1); y++)
//     {
//         t += id[(y + r) * w];
//         t -= id[0];
//         od[y * w] = t * scale;
//     }

//     // main loop
//     for (int y = (r + 1); y < (h - r); y++)
//     {
//         t += id[(y + r) * w];
//         t -= id[((y - r) * w) - w];
//         od[y * w] = t * scale;
//     }

//     // do right edge
//     for (int y = h - r; y < h; y++)
//     {
//         t += id[(h-1) * w];
//         t -= id[((y - r) * w) - w];
//         od[y * w] = t * scale;
//     }
// }

// __global__ void
// d_boxfilter_x_global(float *id, float *od, int w, int h, int r)
// {
//     unsigned int y = blockIdx.x*blockDim.x + threadIdx.x;
//     d_boxfilter_x(&id[y * w], &od[y * w], w, h, r);
// }

// __global__ void
// d_boxfilter_y_global(float *id, float *od, int w, int h, int r)
// {
//     unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
//     d_boxfilter_y(&id[x], &od[x], w, h, r);
// }


double convolutionCUDA(float *output, float *input, float * filter, int type) {
    float *d_input; float *d_output; float * d_output2;float *d_filter;
    cudaError_t err;
    timer kernelTime = timer("kernelTime");
    timer memoryTime = timer("memoryTime");
    int input_size = 0;
    int filter_size = 0;

    if(type > 0){
        input_size = input_height_sep*input_width_sep*sizeof(float);
        filter_size = filter_width_sep * sizeof(float);
    }else{
        input_size = input_height*input_width* sizeof(float);
        filter_size = filter_height*filter_width*sizeof(float);
    }

    // memory allocation

    err = cudaMalloc((void **)&d_input, input_size);
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString( err )); }

    err = cudaMalloc((void **)&d_output2, input_size);
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output2: %s\n",    cudaGetErrorString( err )); }

    err = cudaMalloc((void **)&d_output, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString( err )); }

    err = cudaMalloc((void **)&d_filter, filter_size);


    memoryTime.start();
    // host to device

    err = cudaMemcpy(d_input, input, input_size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device input: %s\n", cudaGetErrorString( err ));  }

    // zero the result array
    err = cudaMemset(d_output, 0, image_height*image_width*sizeof(float));
    err = cudaMemset(d_output2, 0, image_height*image_width*sizeof(float));
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemset output: %s\n", cudaGetErrorString( err ));  }

    if(type == -1){
         err = cudaMemcpyToSymbol(d_kernel, filter, filter_size);
        if (err != cudaSuccess) { fprintf(stderr, "Error in cudaSimbolKernel output: %s\n", cudaGetErrorString( err ));  }
    }
    if(type == 0){
        err = cudaMemcpy(d_filter, filter, filter_size, cudaMemcpyHostToDevice);
        if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device filter: %s\n", cudaGetErrorString( err ));  }
    }



    memoryTime.stop();

    //measure the GPU function
    kernelTime.start();

    if(type == 0){
        dim3 threads(block_size_x, block_size_y);
        dim3 grid(int(ceilf(image_width/(float)threads.x)),
                    int(ceilf(image_height/(float)threads.y)) );
        convolution_kernel_naive<<<grid, threads>>>(d_output, d_input,d_filter);
    }else if(type == -1){
        dim3 threads(block_size_x1, block_size_y1);
        dim3 grid(int(ceilf(image_width/(float)threads.x)),
                    int(ceilf(image_height/(float)threads.y)) );
        convolution_kernel_shared<<<grid, threads>>>(d_output, d_input);
    }else if(type == -2){
        dim3 threads(5, 6);
        dim3 grid(int(ceilf(image_width/(float)threads.x)),
                    int(ceilf(image_height/(float)threads.y)) );
        convolution_kernel_shared<<<grid, threads>>>(d_output, d_input);
    }else if(type == -3){
        dim3 threads(8, 8);
        dim3 grid(int(ceilf(image_width/(float)threads.x)),
                    int(ceilf(image_height/(float)threads.y)) );
        convolution_kernel_shared<<<grid, threads>>>(d_output, d_input);
    }else if(type > 0){
        // dim3 threads_h(block_size_hor,1);
        // dim3 grid_h( int(ceilf(image_width/(float)threads_h.x)), int(image_height));
        // convolution_kernel_horizontal<<<grid_h, threads_h>>>(d_output2, d_input);

        // dim3 threads_v(1,block_size_ver);
        // dim3 grid_v(int(image_width),int(ceilf(image_height/(float)threads_v.y)));
        // convolution_kernel_vertical<<<grid_v, threads_v>>>(d_output, d_output2);
    }
    cudaDeviceSynchronize();
    kernelTime.stop();


    //check to see if all went well
    err = cudaGetLastError();
    if (err != cudaSuccess) { fprintf(stderr, "Error during kernel launch convolution_kernel: %s\n", cudaGetErrorString( err )); }

    //copy the result back to host memory
    memoryTime.start();
    err = cudaMemcpy(output, d_output, image_height*image_width*sizeof(float), cudaMemcpyDeviceToHost);
    memoryTime.stop();
    if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy device to host output: %s\n", cudaGetErrorString( err )); }

    err = cudaFree(d_input);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_input: %s\n", cudaGetErrorString( err )); }
    err = cudaFree(d_output);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString( err )); }
     err = cudaFree(d_output2);
    if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString( err )); }

    if(type == 0){
        err = cudaFree(d_filter);
        if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_filter: %s\n", cudaGetErrorString( err )); }
    }

#ifdef DEBUG
    cout << "convolution (kernel): \t\t" << kernelTime << endl;
    cout << "convolution (memory): \t\t" << memoryTime << endl;
#endif
    return kernelTime.getTimeInSeconds() + memoryTime.getTimeInSeconds();
}

int compare_arrays(float *a1, float *a2, int n) {
    int errors = 0;
    int print = 0;

    for (int i=0; i<n; i++) {

        if (isnan(a1[i]) || isnan(a2[i])) {
            errors++;
            if (print < 10) {
                print++;
                fprintf(stderr, "Error NaN detected at i=%d,\t a1= %10.7e \t a2= \t %10.7e\n",i,a1[i],a2[i]);
            }
        }

        float diff = (a1[i]-a2[i])/a1[i];
        if (diff > 1e-6f) {
            errors++;
            if (print < 10) {
                print++;
                fprintf(stderr, "Error detected at i=%d, \t a1= \t %10.7e \t a2= \t %10.7e \t rel_error=\t %10.7e\n",i,a1[i],a2[i],diff);
            }
        }

    }

    return errors;
}

void testArrays(float *a1, float *a2)
{
    int errors=0;
    errors += compare_arrays(a1, a2, image_height*image_width);
    if (errors > 0) {
        printf("TEST FAILED! %d errors!\n", errors);
    } else {
        printf("TEST PASSED!\n");
    }
}

void runExperiment(float *seq,float *output, float *input, float *filter, int type){
    size_t iterations = 100;
    double time_sum = 0;
    for(size_t i = 0 ; i < iterations;++i){
        time_sum += convolutionCUDA(output, input, filter,type);
    }
    printf("time: %e \n'",time_sum / iterations);
#ifdef DEBUG
    testArrays(seq,output);
#endif
}
int main() {
    int i,j;
    cudaError_t err;

    //allocate arrays and fill them
    float *input = (float *) malloc(input_height * input_width * sizeof(float));
    float *input_sep = (float *) malloc(input_height_sep * input_width_sep * sizeof(float));
    float *output1 = (float *) calloc(image_height * image_width, sizeof(float));
    float *output2 = (float *) calloc(image_height * image_width, sizeof(float));
    float *output3 = (float *) calloc(image_height * image_width, sizeof(float));
    float *output4 = (float *) calloc(image_height * image_width, sizeof(float));

    float * filter = (float *) calloc(filter_height * filter_width , sizeof(float));;


    for (i=0; i< input_height * input_width; i++) {
        input[i] = (float) (i % SEED);
    }
    for(i = 0; i < input_height;++i){
        for(j=0;j<input_width;++j){
            if(i> 0 && i < image_height+1 && j> 0 && j< image_width+1 ){
                 input_sep[(i-1)*(input_width_sep)+j-1] = input[i*input_width + j];
            }
        }
    }


//THis is specific for a W==H smoothening filteri, where W and H are odd.
    for (i=0; i<filter_height * filter_width; i++) {
      filter[i] = 1.0;
    }

    for (i=filter_width+1; i<(filter_height - 1) * filter_width; i++) {
        if (i % filter_width > 0 && i % filter_width < filter_width-1)
            filter[i]+=1.0;
    }

    filter[filter_width*filter_height/2]=3.0;

    for(i = 0; i < filter_height * filter_width; ++i){
        if(i % filter_width == 0)
            printf("\n");
        printf("%e ",filter[i]);

    }
    printf("\n");




//end initialization

    //measure the CPU function
    convolutionSeq(output1, input, filter);
    //measure the GPU function
    printf("DUMMY:\n");
   runExperiment(output1,output2, input, filter,0);

   printf("5x5:\n");
   runExperiment(output1,output2, input, filter,-1);


//    printf("5x6:\n");
//    runExperiment(output1,output2, input, filter,-2);

//    printf("8x8:\n");
//    runExperiment(output1,output2, input, filter,-3);





    //free(input);
    //free(input_sep);
    free(filter);
    free(output1);
    free(output2);
    free(output3);
    free(output4);

    return 0;
}



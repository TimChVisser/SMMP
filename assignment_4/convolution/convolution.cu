#include <stdio.h>
#include <string.h>
#include "timer.h"

// #define image_height 8192
#define image_height 256
#define image_width 8192
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

#define block_size_hor_x 64
#define block_size_hor_y 16

#define block_size_ver_x 1024
#define block_size_ver_y 1
#define DEBUG
#define SEED 1234


#define block_size_x1 16
#define block_size_y1 16



using namespace std;

__constant__ float d_kernel[filter_width*filter_height];


static void checkCudaCall(cudaError_t result) {
    if (result != cudaSuccess) {
        printf("cuda error %s \n",cudaGetErrorString( result ));
        exit(1);
    }
}
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

    if(x>input_width-1 || y > input_height -1){
        return;
    }
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
    if(x<image_width && y < image_height){
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
    }
}




__global__ void
boxfilter_horiz(float * input, float *output, int filter_size, int radius)
{
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int x = blockIdx.x * blockDim.x + threadIdx.x;

    //printf("%i \n",blockDim.x);
    int offset = 0;
    if (blockIdx.x > 0){
        offset = -radius;
    }
    int start = y*input_width + x + offset;
    float t = 0;
    for(int i = 0; i < filter_size; ++ i){
        t+=input[start+i];
    }
    output[start+radius] = t;


    for (int i = 1; i < blockDim.x ; i++)
    {
        t += input[start + i +filter_size -1];
        t -= input[start + i - 1];
        output[start + i+radius] = t;
    }

}

__global__ void
boxfilter_vert(float * input, float *output, int filter_size, int radius)
{
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int x = blockIdx.x * blockDim.x + threadIdx.x;

    // int row = y*width;
    float t = 0;
    int offset = 0;
    if (blockIdx.y > 0){
        offset = -radius;
    }
    for(int i = 0; i < filter_size; ++ i){
        t+=input[(y+i+ offset)*input_width + x];
    }
    output[(offset + radius)* input_width +x] = t;

    for (int i = 1; i < blockDim.y ; i++)
    {
        t += input[(y+ i + filter_size + offset -1)* input_width + x];
        t -= input[(y+ i -1 + offset)* input_width + x];
        output[(y + i + radius+offset)* input_width+x] = t;
    }

}

__global__ void
merge(float * input,float * input1, float * input2, float *output,int radius){
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    size_t idx = (y+ radius) * input_width + x+radius;
    size_t idx_o = y * image_width + x;
    if(y > radius && x > radius && y < image_height+radius && x < image_width + radius)
        output[idx_o]= (input[idx] + input2[idx] + input1[idx]) / 35;

}

// __global__ void convolve(float * input,float * input1, float * input2, float *output,int radius)



double convolutionCUDA(float *output, float *input, float * filter, int type) {
    float *d_input; float *d_output; float * temp; float * res5;float * res3;float *d_filter;

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

    checkCudaCall(cudaMalloc((void **)&d_input, input_size));

    checkCudaCall(cudaMalloc((void **)&temp, input_size));

    checkCudaCall(cudaMalloc((void **)&res5, input_size));

    checkCudaCall(cudaMalloc((void **)&res3, input_size));

    checkCudaCall(cudaMalloc((void **)&d_output, image_height*image_width*sizeof(float)));
    checkCudaCall(cudaMalloc((void **)&d_filter, filter_size));


    memoryTime.start();
    // host to device

    checkCudaCall(cudaMemcpy(d_input, input, input_size, cudaMemcpyHostToDevice));

    // zero the result array
    checkCudaCall(cudaMemset(d_output, 0, image_height*image_width*sizeof(float)));

    if(type == -1){
         checkCudaCall(cudaMemcpyToSymbol(d_kernel, filter, filter_size));

    }
    if(type == 0){
        checkCudaCall(cudaMemcpy(d_filter, filter, filter_size, cudaMemcpyHostToDevice));
    }

    if(type > 0){
        checkCudaCall(cudaMemset(temp, 0, input_size));
        checkCudaCall(cudaMemset(res3, 0, input_size));
        checkCudaCall(cudaMemset(res5, 0, input_size));
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
    }else if(type > 0){
        dim3 threads_h(block_size_hor_x,block_size_hor_y);
        dim3 grid_h( int(ceilf(image_width/(float)threads_h.x)),
                    int(ceilf(image_height/(float)threads_h.y)));


        dim3 threads_v(block_size_ver_x,block_size_ver_y);
        dim3 grid_v(int(ceilf(image_width/(float)threads_v.x)),
                    int(ceilf(image_height/(float)threads_v.y)));

        boxfilter_horiz<<<grid_h, threads_h>>>(temp, d_input,5,2);
        boxfilter_vert<<<grid_v, threads_v>>>(res5, temp,5,2);

        // apply 3x3 box filter
        boxfilter_horiz<<<grid_h, threads_h>>>(temp, d_input,3,1);
        boxfilter_vert<<<grid_v, threads_v>>>(res3, temp,3,1);

        // merge result
        threads_v= dim3(32,32);
        grid_v = dim3(int(ceilf(input_width/(float)threads_v.x)),
                        int(ceilf(input_height/(float)threads_v.y)));
        merge<<<grid_v, threads_v>>>(d_input,res5,res3,d_output, 2);
    }
    cudaDeviceSynchronize();
    kernelTime.stop();


    //check to see if all went well
    checkCudaCall(cudaGetLastError());

    //copy the result back to host memory
    memoryTime.start();
    checkCudaCall(cudaMemcpy(output, d_output, image_height*image_width*sizeof(float),cudaMemcpyDeviceToHost));
    memoryTime.stop();


    checkCudaCall(cudaFree(d_input));
    checkCudaCall(cudaFree(d_output));
     checkCudaCall(cudaFree(temp));
    checkCudaCall(cudaFree(res5));
    checkCudaCall(cudaFree(res3));
    checkCudaCall(cudaFree(d_filter));


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
    size_t iterations = 10;
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

   printf("Shared:\n");
   runExperiment(output1,output3, input, filter,-1);


//    printf("separated:\n");
//     runExperiment(output1,output4, input, filter,1);






    free(input);
    free(input_sep);
    free(filter);
    free(output1);
    free(output2);
    free(output3);
    free(output4);

    return 0;
}



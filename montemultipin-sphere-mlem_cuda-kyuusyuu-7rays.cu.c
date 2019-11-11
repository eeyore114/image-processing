/*3次元シングルピンホール-mlem法*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct{
	int image_size , detecter_width , detecter_hight , detecter_num;
	float radius_rotate , radius_sphere , distance , detecter_pixel , region_pixel , open_degree1 , open_degree2 , open_degree3 , colimater_size;
}general;

typedef struct{
    float x1 , x2 , x3 , y1 ,y2 , y3 , z1 , z2 , z3 , s , t1 , t2 , t3 , u , degree , theta , housen_s , housen_t , housen_u , housen_x , housen_y , housen_z;
} colimater;

__global__ void syokika(float *g_gpu , general *gene_gpu);
__global__ void syokika2(float *h_gpu , general *gene_gpu);
__device__ void plot(float x , float y , float z , general *gene_gpu , int i , int j , int k , float *f_gpu , float *g_gpu , float *MU_gpu , float *mu_l , int *pin_num);
__device__ void kakuritu_plot(float x , float y , float z , general *gene_gpu , int i , int j , int k , float *F_gpu , float *C_gpu , int *pin_num);
__device__ void replot(float x , float y , float z , general *gene_gpu , int i , int j , int k , float *R_gpu , float *h_gpu , int *pin_num);
__global__ void kakuritu(float *F_gpu , float *C_gpu , general *gene_gpu , colimater *coli_data_gpu);
__global__ void mlemprojection(float *f_gpu , float *g_gpu , general *gene_gpu , colimater *coli_data_gpu , float *MU_gpu);
__global__ void rate(float *g_gpu , float *Z_gpu , float *R_gpu , general *gene_gpu);
__global__ void mlemreprojection(float *R_gpu , float *h_gpu , general *gene_gpu , colimater *coli_data_gpu);
__global__ void divide(float *h_gpu , float *f_gpu , float *C_gpu , general *gene_gpu);

int main(){
	int k = 30;//ML-EM法を回す回数
	//general gene = {128 , 512 , 256 , 120 , 25.0f , 10.0f , 7.5f , 0.08f , 0.2f , 45.0f , 15.0f , 30.0f};
	general gene = {128 , 600 , 300 , 120 , 25.0f , 10.0f , 7.5f , 0.08f , 0.2f , 45.0f , 3.0f , 24.0f , 0.25f};
    //general gene = {128 , 600 , 300 , 12 , 25.0f , 10.0f , 7.5f , 0.08f , 0.2f , 39.0f , 9.0f , 24.0f};
	general* gene_gpu;
	float* Z = (float*)calloc(gene.detecter_hight * gene.detecter_width * gene.detecter_num , sizeof(float));//用意された投影データ
	float* f = (float*)calloc(gene.image_size * gene.image_size * gene.image_size , sizeof(float));//全画素1の画像
	float* g = (float*)calloc(gene.detecter_hight * gene.detecter_width * gene.detecter_num , sizeof(float));//全画素1の投影データ
	float* C = (float*)calloc(gene.image_size * gene.image_size * gene.image_size , sizeof(float));//検出確率 
	float* R = (float*)calloc(gene.detecter_hight * gene.detecter_width * gene.detecter_num , sizeof(float));//比
	float* h = (float*)calloc(gene.image_size * gene.image_size * gene.image_size , sizeof(float));
	float* F = (float*)calloc(gene.detecter_hight * gene.detecter_width * gene.detecter_num , sizeof(float));//投影データの大きさ分画素値を1にしたやつ
	float* MU = (float*)malloc(sizeof(float) * gene.image_size * gene.image_size * gene.image_size);//muマップ
    float *Z_gpu , *f_gpu , *g_gpu , *C_gpu , *R_gpu , *h_gpu , *F_gpu , *MU_gpu;
	
    colimater* coli_data = (colimater*)malloc(gene.detecter_num * sizeof(colimater));
    for (int i = 0; i < gene.detecter_num; i++){
        coli_data[i].s = gene.radius_rotate;
        /*coli_data[i].t1 = -13.0f;
        coli_data[i].t2 = 0.0f;
        coli_data[i].t3 = 13.0f;*/
        coli_data[i].t1 = -10.0f;
        coli_data[i].t2 = 0.0f;
        coli_data[i].t3 = 10.0f;
        coli_data[i].u = 0.0f;
        coli_data[i].degree = (360.0 / gene.detecter_num) * i;
        coli_data[i].theta = M_PI * coli_data[i].degree / 180.0f;
        coli_data[i].x1 = coli_data[i].s * cos(coli_data[i].theta) - coli_data[i].t1 * sin(coli_data[i].theta);
        coli_data[i].x2 = coli_data[i].s * cos(coli_data[i].theta) - coli_data[i].t2 * sin(coli_data[i].theta);
        coli_data[i].x3 = coli_data[i].s * cos(coli_data[i].theta) - coli_data[i].t3 * sin(coli_data[i].theta);
        coli_data[i].y1 = coli_data[i].s * sin(coli_data[i].theta) + coli_data[i].t1 * cos(coli_data[i].theta);
        coli_data[i].y2 = coli_data[i].s * sin(coli_data[i].theta) + coli_data[i].t2 * cos(coli_data[i].theta);
        coli_data[i].y3 = coli_data[i].s * sin(coli_data[i].theta) + coli_data[i].t3 * cos(coli_data[i].theta);
        coli_data[i].z1 = coli_data[i].u;
        coli_data[i].z2 = coli_data[i].u;
        coli_data[i].z3 = coli_data[i].u;
        coli_data[i].housen_s = 1.0f;
        coli_data[i].housen_t = 0.0f;
        coli_data[i].housen_u = 0.0f;
        coli_data[i].housen_x = coli_data[i].housen_s * cos(coli_data[i].theta) - coli_data[i].housen_t * sin(coli_data[i].theta);
        coli_data[i].housen_y = coli_data[i].housen_s * sin(coli_data[i].theta) + coli_data[i].housen_t * cos(coli_data[i].theta);
        coli_data[i].housen_z = coli_data[i].housen_u;
    }
    colimater* coli_data_gpu;

    //画像入力
    {
        char fname1[1000] = {};
	    sprintf(fname1 , "multi-sensitivity-newprojectiondata_%d-%d-%d.raw" , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
        //sprintf(fname1 , "multi-sensitivity-newprojectiondata2_%d-%d-%d.raw" , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
        //sprintf(fname1 , "new-sphere-monte-multi_cuda-float_0_%d-%d-%d.raw" , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
	    //sprintf(fname1 , "new-sphere-monte-multi_cuda-float_0_%d-%d-%d.raw" , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
	    //sprintf(fname1 , "multipin-3D-sphere-touei-new_float-kakuninn_%d-%d-%d.raw" , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
	    //sprintf(fname1 , "multipin-3D-sphere-touei_float_%d-%d-%d.raw" , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
	    FILE* fp1 = fopen(fname1 , "rb");
	    fread( Z , sizeof(float) , gene.detecter_hight * gene.detecter_width * gene.detecter_num, fp1 );
	    fclose( fp1 );
    }

    //mumap入力
    {
        FILE* fp2 = fopen( "mumap_H2Ofor140keV_128-128-128.raw" , "rb" );
        fread( MU , sizeof(float) , gene.image_size * gene.image_size * gene.image_size , fp2 );
        fclose( fp2 );
    }
    
    for(int i = 0; i < gene.image_size * gene.image_size * gene.image_size; i++){
    	f[i] = 1.0f;
    }

    for(int i = 0; i < gene.detecter_hight * gene.detecter_width * gene.detecter_num; i++){
    	F[i] = 1.0f;
    }

    // デバイス(GPU)のメモリ領域確保
    cudaMalloc((void**)&Z_gpu , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num);
    cudaMalloc((void**)&f_gpu , sizeof(float) * gene.image_size * gene.image_size * gene.image_size);
    cudaMalloc((void**)&g_gpu , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num);
    cudaMalloc((void**)&C_gpu , sizeof(float) * gene.image_size * gene.image_size * gene.image_size);
    cudaMalloc((void**)&R_gpu , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num);
    cudaMalloc((void**)&h_gpu , sizeof(float) * gene.image_size * gene.image_size * gene.image_size);
    cudaMalloc((void**)&coli_data_gpu , sizeof(colimater) * gene.detecter_num);
    cudaMalloc((void**)&gene_gpu , sizeof(general) * 1);
    cudaMalloc((void**)&F_gpu , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num);
    cudaMalloc((void**)&MU_gpu , sizeof(float) * gene.image_size * gene.image_size * gene.image_size);
    // ホスト(CPU)からデバイス(GPU)へ転送
    cudaMemcpy(Z_gpu , Z , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num, cudaMemcpyHostToDevice);
    cudaMemcpy(f_gpu , f , sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyHostToDevice);
    cudaMemcpy(g_gpu , g , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num , cudaMemcpyHostToDevice);
    cudaMemcpy(C_gpu , C , sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyHostToDevice);
    cudaMemcpy(R_gpu , R , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num , cudaMemcpyHostToDevice);
    cudaMemcpy(h_gpu , h , sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyHostToDevice);
    cudaMemcpy(coli_data_gpu , coli_data , sizeof(colimater) * gene.detecter_num, cudaMemcpyHostToDevice);
    cudaMemcpy(gene_gpu , &gene , sizeof(general) * 1, cudaMemcpyHostToDevice);
    cudaMemcpy(F_gpu , F , sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num , cudaMemcpyHostToDevice);
    cudaMemcpy(MU_gpu, MU , sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyHostToDevice);

    int image_w = gene.detecter_width;
    int image_h = gene.detecter_hight;
    int image_d = gene.detecter_num;
    int thread_num = 2;
    const int blocks_x  = image_w / thread_num;
	const int blocks_y  = image_h / thread_num;
	const int blocks_z  = image_d / thread_num;
	const int threads_x = thread_num;
	const int threads_y = thread_num;
	const int threads_z = thread_num;
	dim3 blocks1 ( blocks_x,  blocks_y, blocks_z);
  	dim3 threads1(threads_x, threads_y, threads_z);
	dim3 blocks(blocks_x,  blocks_y);
	dim3 threads(threads_x, threads_y);
	dim3 blocks2(64 , 64 , 64);
	dim3 threads2(2 , 2 , 2);
	/*dim3 blocks(256 , 128);
	dim3 threads(2 , 2);
	dim3 blocks1(256 , 128 , 60);
	dim3 threads1(2 , 2 , 2);
	dim3 blocks2(64 , 64 , 64);
	dim3 threads2(2 , 2 , 2);*/


	kakuritu <<< blocks , threads >>> (F_gpu , C_gpu , gene_gpu , coli_data_gpu);
	cudaThreadSynchronize();
    // デバイス(GPU)からホスト(CPU)へ転送
    cudaMemcpy(C, C_gpu, sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyDeviceToHost);
    
	for(int kaisuu = 0; kaisuu < k; kaisuu++){
		syokika <<< blocks1 , threads1 >>> (g_gpu , gene_gpu);
		cudaThreadSynchronize();
        // デバイス(GPU)からホスト(CPU)へ転送
        cudaMemcpy(g, g_gpu, sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num, cudaMemcpyDeviceToHost);
        
        syokika2 <<< blocks2 , threads2 >>> (h_gpu , gene_gpu);
        cudaThreadSynchronize();
        // デバイス(GPU)からホスト(CPU)へ転送
        cudaMemcpy(h, h_gpu, sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyDeviceToHost);
        
    	mlemprojection <<< blocks , threads >>> (f_gpu , g_gpu , gene_gpu , coli_data_gpu , MU_gpu);
    	cudaThreadSynchronize();
        // デバイス(GPU)からホスト(CPU)へ転送
        cudaMemcpy(g, g_gpu, sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num, cudaMemcpyDeviceToHost);
        
        rate <<< blocks1 , threads1 >>> (g_gpu , Z_gpu , R_gpu , gene_gpu);
        cudaThreadSynchronize();
        // デバイス(GPU)からホスト(CPU)へ転送
        cudaMemcpy(R, R_gpu, sizeof(float) * gene.detecter_hight * gene.detecter_width * gene.detecter_num, cudaMemcpyDeviceToHost);
        
        mlemreprojection <<< blocks , threads >>> (R_gpu , h_gpu , gene_gpu , coli_data_gpu);
        cudaThreadSynchronize();
        // デバイス(GPU)からホスト(CPU)へ転送
        cudaMemcpy(h, h_gpu, sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyDeviceToHost);
        

        divide <<< blocks2 , threads2 >>> (h_gpu , f_gpu , C_gpu , gene_gpu);
        cudaThreadSynchronize();
        // デバイス(GPU)からホスト(CPU)へ転送
        cudaMemcpy(f, f_gpu, sizeof(float) * gene.image_size * gene.image_size * gene.image_size, cudaMemcpyDeviceToHost);
    }

    //画像出力
 	{
    	char fname1[1000] = {};
    	//sprintf(fname1, "otameshimlem-%d_%d-%d-%d.raw", k , gene.detecter_width , gene.detecter_hight , gene.detecter_num);
    	//sprintf(fname1, "otameshimlem-%d_%d-%d-%d.raw", k , gene.image_size , gene.image_size , gene.image_size);
    	//sprintf(fname1, "multipin-3D-shepp-mlem_float1-%d_%d-%d-%d.raw", k , gene.image_size , gene.image_size , gene.image_size);
    	sprintf(fname1, "multipin-3D-sphere-kyuusyuu-kanndo-7rays-monte-120touei_float-%d_%d-%d-%d.raw", k , gene.image_size , gene.image_size , gene.image_size);
    	FILE* fp3 = fopen(fname1, "wb");
		fwrite( f, sizeof(float) , gene.image_size * gene.image_size * gene.image_size , fp3);
		//fwrite( g, sizeof(float) , gene.detecter_width * gene.detecter_hight * gene.detecter_num , fp3);
 		fclose(fp3);
 	}
	
 	free(Z);
 	free(f);
 	free(g);
 	free(C);
 	free(R);
 	free(h);
 	free(F);
 	cudaFree(Z_gpu);
 	cudaFree(f_gpu);
 	cudaFree(g_gpu);
 	cudaFree(C_gpu);
 	cudaFree(R_gpu);
 	cudaFree(h_gpu);
 	cudaFree(F_gpu);
 	return 0;
}

__global__ void syokika(float *g_gpu , general *gene_gpu){
	general gene_gpu1 = gene_gpu[0];
	int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    g_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] = 0.0f;
}

__global__ void syokika2(float *h_gpu , general *gene_gpu){
	general gene_gpu1 = gene_gpu[0];
	int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    h_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k] = 0.0f;
}

__device__ void kakuritu_plot(float x , float y , float z , general *gene_gpu , int i , int j , int k , float *F_gpu , float *C_gpu , int *pin_num){
	general gene_gpu1 = gene_gpu[0];
	float X = x / gene_gpu1.region_pixel;
	float Y = y / gene_gpu1.region_pixel;
	float Z = z / gene_gpu1.region_pixel;
		    			
	if(-(gene_gpu1.image_size - 1.0f) / 2.0f < X && X < (gene_gpu1.image_size - 1.0f) / 2.0f && -(gene_gpu1.image_size - 1.0f) / 2.0f < Y && Y < (gene_gpu1.image_size - 1.0f)/2.0f && -(gene_gpu1.image_size - 1.0f) / 2.0f < Z && Z < (gene_gpu1.image_size - 1.0f)/2.0f){
		float x0 = ceil(X - 0.5f) + 0.5f;
		float y0 = floor(Y - 0.5f) + 0.5f;
		float z0 = ceil(Z - 0.5f) + 0.5f;
									      
		int I = (gene_gpu1.image_size - 1.0f) / 2.0f + x0;
		int J = (gene_gpu1.image_size - 1.0f) / 2.0f + y0;
		int K = (gene_gpu1.image_size - 1.0f) / 2.0f - z0;
						      
		float dx = fabs(X - x0);
		float dy = fabs(Y - y0);
		float dz = fabs(Z - z0);
									      
		float s1 = dx * dy * dz;
		float s2 = dx * (1.0f - dy) * dz;
		float s3 = (1.0f - dx) * dy * dz;
		float s4 = (1.0f - dx) * (1.0f - dy) * dz;
		float s5 = dx * dy * (1.0f - dz);
		float s6 = dx * (1.0f - dy) * (1.0f - dz);
		float s7 = (1.0f - dx) * dy * (1.0f - dz);
		float s8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);
								      
		int f1 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * I + J;
		int f2 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * I + (J + 1);
		int f3 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * (I - 1) + J;
		int f4 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * (I - 1) + (J + 1);
		int f5 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * I + J;
		int f6 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * I + (J + 1);
		int f7 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * (I - 1) + J;
		int f8 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * (I - 1) + (J + 1);

        if(*pin_num == 3){
            atomicAdd(&C_gpu[f1] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s8) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f2] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s7) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f3] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s6) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f4] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s5) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f5] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s4) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f6] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s3) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f7] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s2) * (1.0f/4.0f));   
            atomicAdd(&C_gpu[f8] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s1) * (1.0f/4.0f));   
        }
		else{						      
    		atomicAdd(&C_gpu[f1] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s8) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f2] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s7) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f3] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s6) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f4] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s5) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f5] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s4) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f6] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s3) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f7] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s2) * (1.0f/8.0f));	
    		atomicAdd(&C_gpu[f8] , (F_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s1) * (1.0f/8.0f));	
	    }
    }
}

__device__ void plot(float x , float y , float z , general *gene_gpu , int i , int j , int k , float *f_gpu , float *g_gpu , float *MU_gpu , float *mu_l , int *pin_num){
	general gene_gpu1 = gene_gpu[0];
	float val = 0.0f;
    float sum = 0.0f;
    float X = x / gene_gpu1.region_pixel;
	float Y = y / gene_gpu1.region_pixel;
	float Z = z / gene_gpu1.region_pixel;
	if(-(gene_gpu1.image_size - 1.0f) / 2.0f < X && X < (gene_gpu1.image_size - 1.0f) / 2.0f && -(gene_gpu1.image_size - 1.0f) / 2.0f < Y && Y < (gene_gpu1.image_size - 1.0f)/2.0f && -(gene_gpu1.image_size - 1.0f) / 2.0f < Z && Z < (gene_gpu1.image_size - 1.0f)/2.0f){
		float x0 = ceil(X - 0.5f) + 0.5f;
		float y0 = floor(Y - 0.5f) + 0.5f;
		float z0 = ceil(Z - 0.5f) + 0.5f;
									      
		int I = (gene_gpu1.image_size - 1.0f) / 2.0f + x0;
		int J = (gene_gpu1.image_size - 1.0f) / 2.0f + y0;
		int K = (gene_gpu1.image_size - 1.0f) / 2.0f - z0;
						      
		float dx = fabs(X - x0);
		float dy = fabs(Y - y0);
		float dz = fabs(Z - z0);
									      
		float s1 = dx * dy * dz;
		float s2 = dx * (1.0f - dy) * dz;
		float s3 = (1.0f - dx) * dy * dz;
		float s4 = (1.0f - dx) * (1.0f - dy) * dz;
		float s5 = dx * dy * (1.0f - dz);
		float s6 = dx * (1.0f - dy) * (1.0f - dz);
		float s7 = (1.0f - dx) * dy * (1.0f - dz);
		float s8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);
								      
		int f1 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * I + J;
		int f2 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * I + (J + 1);
		int f3 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * (I - 1) + J;
		int f4 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * (I - 1) + (J + 1);
		int f5 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * I + J;
		int f6 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * I + (J + 1);
		int f7 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * (I - 1) + J;
		int f8 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * (I - 1) + (J + 1);

        if(*pin_num == 3){
            val = (f_gpu[f1] * s8 + f_gpu[f2] * s7 + f_gpu[f3] * s6 + f_gpu[f4] * s5 + f_gpu[f5] * s4 + f_gpu[f6] * s3 + f_gpu[f7] * s2 + f_gpu[f8] * s1) * (1.0f/4.0f);
        }
        else{
            val = (f_gpu[f1] * s8 + f_gpu[f2] * s7 + f_gpu[f3] * s6 + f_gpu[f4] * s5 + f_gpu[f5] * s4 + f_gpu[f6] * s3 + f_gpu[f7] * s2 + f_gpu[f8] * s1) * (1.0f/8.0f);
        }
		*mu_l += MU_gpu[f1] * s8 + MU_gpu[f2] * s7 + MU_gpu[f3] * s6 + MU_gpu[f4] * s5 + MU_gpu[f5] * s4 + MU_gpu[f6] * s3 + MU_gpu[f7] * s2 + MU_gpu[f8] * s1;
        /*if(i == 91 && j == 130){
            printf("mu_l=%f\n",*mu_l);
        }*/
        sum = val * exp(-*mu_l);
        atomicAdd(&g_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] , sum);	
	}
}

__device__ void replot(float x , float y , float z , general *gene_gpu , int i , int j , int k , float *R_gpu , float *h_gpu , int *pin_num){
	general gene_gpu1 = gene_gpu[0];
	float X = x / gene_gpu1.region_pixel;
	float Y = y / gene_gpu1.region_pixel;
	float Z = z / gene_gpu1.region_pixel;
		    			
	if(-(gene_gpu1.image_size - 1.0f) / 2.0f < X && X < (gene_gpu1.image_size - 1.0f) / 2.0f && -(gene_gpu1.image_size - 1.0f) / 2.0f < Y && Y < (gene_gpu1.image_size - 1.0f)/2.0f && -(gene_gpu1.image_size - 1.0f) / 2.0f < Z && Z < (gene_gpu1.image_size - 1.0f)/2.0f){
		float x0 = ceil(X - 0.5f) + 0.5f;
		float y0 = floor(Y - 0.5f) + 0.5f;
		float z0 = ceil(Z - 0.5f) + 0.5f;
									      
		int I = (gene_gpu1.image_size - 1.0f) / 2.0f + x0;
		int J = (gene_gpu1.image_size - 1.0f) / 2.0f + y0;
		int K = (gene_gpu1.image_size - 1.0f) / 2.0f - z0;
						      
		float dx = fabs(X - x0);
		float dy = fabs(Y - y0);
		float dz = fabs(Z - z0);
									      
		float s1 = dx * dy * dz;
		float s2 = dx * (1.0f - dy) * dz;
		float s3 = (1.0f - dx) * dy * dz;
		float s4 = (1.0f - dx) * (1.0f - dy) * dz;
		float s5 = dx * dy * (1.0f - dz);
		float s6 = dx * (1.0f - dy) * (1.0f - dz);
		float s7 = (1.0f - dx) * dy * (1.0f - dz);
		float s8 = (1.0f - dx) * (1.0f - dy) * (1.0f - dz);
								      
		int f1 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * I + J;
		int f2 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * I + (J + 1);
		int f3 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * (I - 1) + J;
		int f4 = gene_gpu1.image_size * gene_gpu1.image_size * K + gene_gpu1.image_size * (I - 1) + (J + 1);
		int f5 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * I + J;
		int f6 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * I + (J + 1);
		int f7 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * (I - 1) + J;
		int f8 = gene_gpu1.image_size * gene_gpu1.image_size * (K + 1) + gene_gpu1.image_size * (I - 1) + (J + 1);

        // ----------------------------------------------------------
        // float weight = *pin_num == 3 ? 1. / 4. : 1. / 8.;
        // atomicAdd(&h_gpu[f1] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s8) * weight); 
        // atomicAdd(&h_gpu[f2] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s7) * weight); 
        // atomicAdd(&h_gpu[f3] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s6) * weight); 
        // atomicAdd(&h_gpu[f4] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s5) * weight); 
        // atomicAdd(&h_gpu[f5] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s4) * weight); 
        // atomicAdd(&h_gpu[f6] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s3) * weight); 
        // atomicAdd(&h_gpu[f7] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s2) * weight); 
        // atomicAdd(&h_gpu[f8] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s1) * weight);
        // ----------------------------------------------------------
       



		if(*pin_num == 3){
            atomicAdd(&h_gpu[f1] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s8) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f2] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s7) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f3] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s6) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f4] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s5) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f5] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s4) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f6] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s3) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f7] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s2) * (1.0f/4.0f)); 
            atomicAdd(&h_gpu[f8] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s1) * (1.0f/4.0f));
        }
        else{						      
    		atomicAdd(&h_gpu[f1] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s8) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f2] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s7) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f3] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s6) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f4] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s5) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f5] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s4) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f6] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s3) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f7] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s2) * (1.0f/8.0f));	
    		atomicAdd(&h_gpu[f8] , (R_gpu[j + gene_gpu1.detecter_width * i + gene_gpu1.detecter_hight * gene_gpu1.detecter_width * k] * s1) * (1.0f/8.0f));
        }	
	}
}

__global__ void kakuritu(float *F_gpu , float *C_gpu , general *gene_gpu , colimater *coli_data_gpu){
	general gene_gpu1 = gene_gpu[0];
    int pin_num = 0;
    for(int k = 0; k < gene_gpu1.detecter_num; k++){
		int j = blockIdx.x * blockDim.x + threadIdx.x;//横
	    int i = blockIdx.y * blockDim.y + threadIdx.y;//縦
	    /*検出器の座標(s,t,u座標)*/
		float s = gene_gpu1.radius_rotate + gene_gpu1.distance;
		float t = -(((gene_gpu1.detecter_width - 1.0f) / 2.0f) * gene_gpu1.detecter_pixel) + j * gene_gpu1.detecter_pixel;
		float u = (((gene_gpu1.detecter_hight - 1.0f) / 2.0f) * gene_gpu1.detecter_pixel) - i * gene_gpu1.detecter_pixel;
		float detecter_x = s * cos(coli_data_gpu[k].theta) - t * sin(coli_data_gpu[k].theta);
		float detecter_y = s * sin(coli_data_gpu[k].theta) + t * cos(coli_data_gpu[k].theta);
		float detecter_z = u;
		/*どのコリメータの条件に適しているのかの判定*/
		if(t < coli_data_gpu[k].t1){
            for(pin_num = 0; pin_num <= 2; pin_num++){
                if(pin_num == 0 || pin_num == 2){
                    coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                    coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                    if(pin_num == 0){
                        coli_data_gpu[i].z1 = coli_data_gpu[i].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                    }
                    else{
                        coli_data_gpu[i].z1 = coli_data_gpu[i].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                    }
                }
                else{
                    coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                    coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                    coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                }
    			float vector_x1 = coli_data_gpu[k].x1 - detecter_x;
    			float vector_y1 = coli_data_gpu[k].y1 - detecter_y;
    			float vector_z1 = coli_data_gpu[k].z1 - detecter_z;
    			float numerator1 = (-vector_x1) * coli_data_gpu[k].housen_x + (-vector_y1) * coli_data_gpu[k].housen_y + (-vector_z1) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
                float denominator1 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x1 , 2.0f) + pow(-vector_y1 , 2.0f) + pow(-vector_z1 , 2.0f))));
                double in_theta1 = acos(numerator1 / denominator1);
                if((in_theta1 * 180.0f / M_PI) <= gene_gpu1.open_degree1){
                    for(int l = 0; l < 1000; l++){
    					/*単位ベクトル*/
    					float vector_size = sqrt(pow(vector_x1 , 2.0f) + pow(vector_y1 , 2.0f) + pow(vector_z1 , 2.0f));
    					float unitvector_x = vector_x1 / vector_size;
    					float unitvector_y = vector_y1 / vector_size;
    					float unitvector_z = vector_z1 / vector_size;
    					float x = coli_data_gpu[k].x1 + gene_gpu1.region_pixel * l * unitvector_x;
    					float y = coli_data_gpu[k].y1 + gene_gpu1.region_pixel * l * unitvector_y;
    					float z = coli_data_gpu[k].z1 + gene_gpu1.region_pixel * l * unitvector_z;
    					kakuritu_plot(x , y , z , gene_gpu , i , j , k , F_gpu , C_gpu , &pin_num);
    				}
                }
            }
        }
		else{
			if(t < coli_data_gpu[k].t2){
                for(pin_num = 0; pin_num <= 6; pin_num++){
                    if(pin_num == 0 || pin_num == 2){
                        coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        if(pin_num == 0){
                            coli_data_gpu[i].z2 = coli_data_gpu[i].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                        else{
                            coli_data_gpu[i].z2 = coli_data_gpu[i].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                    }
                    else if(pin_num == 1){
                        coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                    }
                    else if(pin_num == 4 || pin_num == 6){
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        if(pin_num == 4){
                            coli_data_gpu[k].z1 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                        else{
                            coli_data_gpu[k].z1 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                    }
                    else if(pin_num == 5){
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                    }
                    else{
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t1 * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t1 * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                    }
    				float vector_x1 = coli_data_gpu[k].x1 - detecter_x;
    				float vector_y1 = coli_data_gpu[k].y1 - detecter_y;
    				float vector_z1 = coli_data_gpu[k].z1 - detecter_z;
    				float numerator1 = (-vector_x1) * coli_data_gpu[k].housen_x + (-vector_y1) * coli_data_gpu[k].housen_y + (-vector_z1) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    	            float denominator1 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x1 , 2.0f) + pow(-vector_y1 , 2.0f) + pow(-vector_z1 , 2.0f))));
    	            float in_theta1 = acos(numerator1 / denominator1);
    	            float vector_x2 = coli_data_gpu[k].x2 - detecter_x;
    				float vector_y2 = coli_data_gpu[k].y2 - detecter_y;
    				float vector_z2 = coli_data_gpu[k].z2 - detecter_z;
    				float numerator2 = (-vector_x2) * coli_data_gpu[k].housen_x + (-vector_y2) * coli_data_gpu[k].housen_y + (-vector_z2) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    	            float denominator2 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x2 , 2.0f) + pow(-vector_y2 , 2.0f) + pow(-vector_z2 , 2.0f))));
    	            float in_theta2 = acos(numerator2 / denominator2);
    	            if((in_theta1 * 180.0f / M_PI) <= gene_gpu1.open_degree2){
    	                for(int l = 0; l < 1000; l++){
    						/*単位ベクトル*/
    						float vector_size = sqrt(pow(vector_x1 , 2.0f) + pow(vector_y1 , 2.0f) + pow(vector_z1 , 2.0f));
    						float unitvector_x = vector_x1 / vector_size;
    						float unitvector_y = vector_y1 / vector_size;
    						float unitvector_z = vector_z1 / vector_size;
    						float x = coli_data_gpu[k].x1 + gene_gpu1.region_pixel * l * unitvector_x;
    						float y = coli_data_gpu[k].y1 + gene_gpu1.region_pixel * l * unitvector_y;
    						float z = coli_data_gpu[k].z1 + gene_gpu1.region_pixel * l * unitvector_z;
    						kakuritu_plot(x , y , z , gene_gpu , i , j , k , F_gpu , C_gpu , &pin_num);
    					}
    	            }
    	            else{
    	                if((in_theta2 * 180.0f / M_PI) <= gene_gpu1.open_degree3){
    	                    for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x2 , 2.0f) + pow(vector_y2 , 2.0f) + pow(vector_z2 , 2.0f));
    							float unitvector_x = vector_x2 / vector_size;
    							float unitvector_y = vector_y2 / vector_size;
    							float unitvector_z = vector_z2 / vector_size;
    							float x = coli_data_gpu[k].x2 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y2 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z2 + gene_gpu1.region_pixel * l * unitvector_z;
    							kakuritu_plot(x , y , z , gene_gpu , i , j , k , F_gpu , C_gpu , &pin_num);
    						}
    	                }
    	            }
                }
	        }
	        else{
	            if(t < coli_data_gpu[k].t3){
                    for(pin_num = 0; pin_num <= 6; pin_num++){
                        if(pin_num == 0 || pin_num == 2){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 0){
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 1){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
                        else if(pin_num == 4 || pin_num == 6){
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 4){
                                coli_data_gpu[k].z2 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z2 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 5){
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                        }
                        else{
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t2 * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t2 * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                        }
    	                float vector_x2 = coli_data_gpu[k].x2 - detecter_x;
    					float vector_y2 = coli_data_gpu[k].y2 - detecter_y;
    					float vector_z2 = coli_data_gpu[k].z2 - detecter_z;
    					float numerator2 = (-vector_x2) * coli_data_gpu[k].housen_x + (-vector_y2) * coli_data_gpu[k].housen_y + (-vector_z2) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator2 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x2 , 2.0f) + pow(-vector_y2 , 2.0f) + pow(-vector_z2 , 2.0f))));
    	                float in_theta2 = acos(numerator2 / denominator2);
    		            float vector_x3 = coli_data_gpu[k].x3 - detecter_x;
    					float vector_y3 = coli_data_gpu[k].y3 - detecter_y;
    					float vector_z3 = coli_data_gpu[k].z3 - detecter_z;
    					float numerator3 = (-vector_x3) * coli_data_gpu[k].housen_x + (-vector_y3) * coli_data_gpu[k].housen_y + (-vector_z3) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator3 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x3 , 2.0f) + pow(-vector_y3 , 2.0f) + pow(-vector_z3 , 2.0f))));
    		            float in_theta3 = acos(numerator3 / denominator3);
    		            if((in_theta2 * 180.0f / M_PI) <= gene_gpu1.open_degree3){
    		                for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x2 , 2.0f) + pow(vector_y2 , 2.0f) + pow(vector_z2 , 2.0f));
    							float unitvector_x = vector_x2 / vector_size;
    							float unitvector_y = vector_y2 / vector_size;
    							float unitvector_z = vector_z2 / vector_size;
    							float x = coli_data_gpu[k].x2 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y2 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z2 + gene_gpu1.region_pixel * l * unitvector_z;
    							kakuritu_plot(x , y , z , gene_gpu , i , j , k , F_gpu , C_gpu , &pin_num);
    						}
    		            }
    		            else{
    		                if((in_theta3 * 180.0f / M_PI) <= gene_gpu1.open_degree2){
    		                    for(int l = 0; l < 1000; l++){
    								/*単位ベクトル*/
    								float vector_size = sqrt(pow(vector_x3 , 2.0f) + pow(vector_y3 , 2.0f) + pow(vector_z3 , 2.0f));
    								float unitvector_x = vector_x3 / vector_size;
    								float unitvector_y = vector_y3 / vector_size;
    								float unitvector_z = vector_z3 / vector_size;
    								float x = coli_data_gpu[k].x3 + gene_gpu1.region_pixel * l * unitvector_x;
    								float y = coli_data_gpu[k].y3 + gene_gpu1.region_pixel * l * unitvector_y;
    								float z = coli_data_gpu[k].z3 + gene_gpu1.region_pixel * l * unitvector_z;
    								kakuritu_plot(x , y , z , gene_gpu , i , j , k , F_gpu , C_gpu , &pin_num);
    							}
    		                }
    		            }
                    }
	            }
	            else{
                    for(pin_num = 3; pin_num <= 6; pin_num++){
                        if(pin_num == 4 || pin_num == 6){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 4){
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 5){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
                        else{
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t3 * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t3 * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
    	                float vector_x3 = coli_data_gpu[k].x3 - detecter_x;
    					float vector_y3 = coli_data_gpu[k].y3 - detecter_y;
    					float vector_z3 = coli_data_gpu[k].z3 - detecter_z;
    					float numerator3 = (-vector_x3) * coli_data_gpu[k].housen_x + (-vector_y3) * coli_data_gpu[k].housen_y + (-vector_z3) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator3 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x3 , 2.0f) + pow(-vector_y3 , 2.0f) + pow(-vector_z3 , 2.0f))));
    		            float in_theta3 = acos(numerator3 / denominator3);
    		            if((in_theta3 * 180.0f / M_PI) <= gene_gpu1.open_degree1){
    		                for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x3 , 2.0f) + pow(vector_y3 , 2.0f) + pow(vector_z3 , 2.0f));
    							float unitvector_x = vector_x3 / vector_size;
    							float unitvector_y = vector_y3 / vector_size;
    							float unitvector_z = vector_z3 / vector_size;
    							float x = coli_data_gpu[k].x3 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y3 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z3 + gene_gpu1.region_pixel * l * unitvector_z;
    							kakuritu_plot(x , y , z , gene_gpu , i , j , k , F_gpu , C_gpu , &pin_num);
    						}
    		            }
                    }
	            }
	        }
	    }
    }
}

__global__ void mlemprojection(float *f_gpu , float *g_gpu , general *gene_gpu , colimater *coli_data_gpu , float *MU_gpu){
	general gene_gpu1 = gene_gpu[0];
    int pin_num = 0;
    for(int k = 0; k < gene_gpu1.detecter_num; k++){
		int j = blockIdx.x * blockDim.x + threadIdx.x;//横
	    int i = blockIdx.y * blockDim.y + threadIdx.y;//縦
	    /*検出器の座標(s,t,u座標)*/
		float s = gene_gpu1.radius_rotate + gene_gpu1.distance;
		float t = -(((gene_gpu1.detecter_width - 1.0f) / 2.0f) * gene_gpu1.detecter_pixel) + j * gene_gpu1.detecter_pixel;
		float u = (((gene_gpu1.detecter_hight - 1.0f) / 2.0f) * gene_gpu1.detecter_pixel) - i * gene_gpu1.detecter_pixel;
		float detecter_x = s * cos(coli_data_gpu[k].theta) - t * sin(coli_data_gpu[k].theta);
		float detecter_y = s * sin(coli_data_gpu[k].theta) + t * cos(coli_data_gpu[k].theta);
		float detecter_z = u;
        float mu_l;
        /*どのコリメータの条件に適しているのかの判定*/
		if(t < coli_data_gpu[k].t1){
            for(pin_num = 0; pin_num <= 2; pin_num++){
                if(pin_num == 0 || pin_num == 2){
                    coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                    coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                    if(pin_num == 0){
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                    }
                    else{
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                    }
                }
                else{
                    coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                    coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                    coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                }

    			float vector_x1 = coli_data_gpu[k].x1 - detecter_x;
    			float vector_y1 = coli_data_gpu[k].y1 - detecter_y;
    			float vector_z1 = coli_data_gpu[k].z1 - detecter_z;
    			float numerator1 = (-vector_x1) * coli_data_gpu[k].housen_x + (-vector_y1) * coli_data_gpu[k].housen_y + (-vector_z1) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
                float denominator1 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x1 , 2.0f) + pow(-vector_y1 , 2.0f) + pow(-vector_z1 , 2.0f))));
                double in_theta1 = acos(numerator1 / denominator1);
                if((in_theta1 * 180.0f / M_PI) <= gene_gpu1.open_degree1){
                    mu_l = 0.0f;
                    for(int l = 0; l < 1000; l++){
    					/*単位ベクトル*/
    					float vector_size = sqrt(pow(vector_x1 , 2.0f) + pow(vector_y1 , 2.0f) + pow(vector_z1 , 2.0f));
    					float unitvector_x = vector_x1 / vector_size;
    					float unitvector_y = vector_y1 / vector_size;
    					float unitvector_z = vector_z1 / vector_size;
    					float x = coli_data_gpu[k].x1 + gene_gpu1.region_pixel * l * unitvector_x;
    					float y = coli_data_gpu[k].y1 + gene_gpu1.region_pixel * l * unitvector_y;
    					float z = coli_data_gpu[k].z1 + gene_gpu1.region_pixel * l * unitvector_z;
    					plot(x , y , z , gene_gpu , i , j , k , f_gpu , g_gpu , MU_gpu , &mu_l , &pin_num);
    				}
                }
            }
        }
		else{
			if(t < coli_data_gpu[k].t2){
                for(pin_num = 0; pin_num <= 6; pin_num++){
                    if(pin_num == 0 || pin_num == 2){
                        coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        if(pin_num == 0){
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                        else{
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                    }
                    else if(pin_num == 1){
                        coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                    }
                    else if(pin_num == 4 || pin_num == 6){
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        if(pin_num == 4){
                            coli_data_gpu[k].z1 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                        else{
                            coli_data_gpu[k].z1 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                    }
                    else if(pin_num == 5){
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                    }
                    else{
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t1 * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t1 * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                    }
    				float vector_x1 = coli_data_gpu[k].x1 - detecter_x;
    				float vector_y1 = coli_data_gpu[k].y1 - detecter_y;
    				float vector_z1 = coli_data_gpu[k].z1 - detecter_z;
    				float numerator1 = (-vector_x1) * coli_data_gpu[k].housen_x + (-vector_y1) * coli_data_gpu[k].housen_y + (-vector_z1) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    	            float denominator1 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x1 , 2.0f) + pow(-vector_y1 , 2.0f) + pow(-vector_z1 , 2.0f))));
    	            float in_theta1 = acos(numerator1 / denominator1);
    	            float vector_x2 = coli_data_gpu[k].x2 - detecter_x;
    				float vector_y2 = coli_data_gpu[k].y2 - detecter_y;
    				float vector_z2 = coli_data_gpu[k].z2 - detecter_z;
    				float numerator2 = (-vector_x2) * coli_data_gpu[k].housen_x + (-vector_y2) * coli_data_gpu[k].housen_y + (-vector_z2) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    	            float denominator2 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x2 , 2.0f) + pow(-vector_y2 , 2.0f) + pow(-vector_z2 , 2.0f))));
    	            float in_theta2 = acos(numerator2 / denominator2);
    	            if((in_theta1 * 180.0f / M_PI) <= gene_gpu1.open_degree2){
                        mu_l = 0.0f;
    	                for(int l = 0; l < 1000; l++){
    						/*単位ベクトル*/
    						float vector_size = sqrt(pow(vector_x1 , 2.0f) + pow(vector_y1 , 2.0f) + pow(vector_z1 , 2.0f));
    						float unitvector_x = vector_x1 / vector_size;
    						float unitvector_y = vector_y1 / vector_size;
    						float unitvector_z = vector_z1 / vector_size;
    						float x = coli_data_gpu[k].x1 + gene_gpu1.region_pixel * l * unitvector_x;
    						float y = coli_data_gpu[k].y1 + gene_gpu1.region_pixel * l * unitvector_y;
    						float z = coli_data_gpu[k].z1 + gene_gpu1.region_pixel * l * unitvector_z;
    						plot(x , y , z , gene_gpu , i , j , k , f_gpu , g_gpu , MU_gpu , &mu_l , &pin_num);
    					}
    	            }
    	            else{
    	                if((in_theta2 * 180.0f / M_PI) <= gene_gpu1.open_degree3){
                            mu_l = 0.0f;
    	                    for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x2 , 2.0f) + pow(vector_y2 , 2.0f) + pow(vector_z2 , 2.0f));
    							float unitvector_x = vector_x2 / vector_size;
    							float unitvector_y = vector_y2 / vector_size;
    							float unitvector_z = vector_z2 / vector_size;
    							float x = coli_data_gpu[k].x2 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y2 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z2 + gene_gpu1.region_pixel * l * unitvector_z;
    							plot(x , y , z , gene_gpu , i , j , k , f_gpu , g_gpu , MU_gpu , &mu_l , &pin_num);
                                /*if(i == 150 && j == 300 && k == 63){
                                    printf("l=%d,mu_l=%f\n",l,&mu_l);
                                }*/
    						}
    	                }
    	            }
                }
	        }
	        else{
	            if(t < coli_data_gpu[k].t3){
                    for(pin_num = 0; pin_num <= 6; pin_num++){
                        if(pin_num == 0 || pin_num == 2){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 0){
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 1){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
                        else if(pin_num == 4 || pin_num == 6){
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 4){
                                coli_data_gpu[k].z2 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z2 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 5){
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                        }
                        else{
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t2 * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t2 * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                        }
    	                float vector_x2 = coli_data_gpu[k].x2 - detecter_x;
    					float vector_y2 = coli_data_gpu[k].y2 - detecter_y;
    					float vector_z2 = coli_data_gpu[k].z2 - detecter_z;
    					float numerator2 = (-vector_x2) * coli_data_gpu[k].housen_x + (-vector_y2) * coli_data_gpu[k].housen_y + (-vector_z2) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator2 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x2 , 2.0f) + pow(-vector_y2 , 2.0f) + pow(-vector_z2 , 2.0f))));
    	                float in_theta2 = acos(numerator2 / denominator2);
    		            float vector_x3 = coli_data_gpu[k].x3 - detecter_x;
    					float vector_y3 = coli_data_gpu[k].y3 - detecter_y;
    					float vector_z3 = coli_data_gpu[k].z3 - detecter_z;
    					float numerator3 = (-vector_x3) * coli_data_gpu[k].housen_x + (-vector_y3) * coli_data_gpu[k].housen_y + (-vector_z3) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator3 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x3 , 2.0f) + pow(-vector_y3 , 2.0f) + pow(-vector_z3 , 2.0f))));
    		            float in_theta3 = acos(numerator3 / denominator3);
    		            if((in_theta2 * 180.0f / M_PI) <= gene_gpu1.open_degree3){
                            mu_l = 0.0f;
    		                for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x2 , 2.0f) + pow(vector_y2 , 2.0f) + pow(vector_z2 , 2.0f));
    							float unitvector_x = vector_x2 / vector_size;
    							float unitvector_y = vector_y2 / vector_size;
    							float unitvector_z = vector_z2 / vector_size;
    							float x = coli_data_gpu[k].x2 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y2 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z2 + gene_gpu1.region_pixel * l * unitvector_z;
    							plot(x , y , z , gene_gpu , i , j , k , f_gpu , g_gpu , MU_gpu , &mu_l , &pin_num);
    						}
    		            }
    		            else{
    		                if((in_theta3 * 180.0f / M_PI) <= gene_gpu1.open_degree2){
                                mu_l = 0.0f;
    		                    for(int l = 0; l < 1000; l++){
    								/*単位ベクトル*/
    								float vector_size = sqrt(pow(vector_x3 , 2.0f) + pow(vector_y3 , 2.0f) + pow(vector_z3 , 2.0f));
    								float unitvector_x = vector_x3 / vector_size;
    								float unitvector_y = vector_y3 / vector_size;
    								float unitvector_z = vector_z3 / vector_size;
    								float x = coli_data_gpu[k].x3 + gene_gpu1.region_pixel * l * unitvector_x;
    								float y = coli_data_gpu[k].y3 + gene_gpu1.region_pixel * l * unitvector_y;
    								float z = coli_data_gpu[k].z3 + gene_gpu1.region_pixel * l * unitvector_z;
    								plot(x , y , z , gene_gpu , i , j , k , f_gpu , g_gpu , MU_gpu , &mu_l , &pin_num);
    							}
    		                }
    		            }
                    }
	            }
	            else{
                    for(pin_num = 3; pin_num <= 6; pin_num++){
                        if(pin_num == 4 || pin_num == 6){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 4){
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 5){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
                        else{
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t3 * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t3 * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
    	                float vector_x3 = coli_data_gpu[k].x3 - detecter_x;
    					float vector_y3 = coli_data_gpu[k].y3 - detecter_y;
    					float vector_z3 = coli_data_gpu[k].z3 - detecter_z;
    					float numerator3 = (-vector_x3) * coli_data_gpu[k].housen_x + (-vector_y3) * coli_data_gpu[k].housen_y + (-vector_z3) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator3 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x3 , 2.0f) + pow(-vector_y3 , 2.0f) + pow(-vector_z3 , 2.0f))));
    		            float in_theta3 = acos(numerator3 / denominator3);
    		            if((in_theta3 * 180.0f / M_PI) <= gene_gpu1.open_degree1){
                            mu_l = 0.0f;
    		                for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x3 , 2.0f) + pow(vector_y3 , 2.0f) + pow(vector_z3 , 2.0f));
    							float unitvector_x = vector_x3 / vector_size;
    							float unitvector_y = vector_y3 / vector_size;
    							float unitvector_z = vector_z3 / vector_size;
    							float x = coli_data_gpu[k].x3 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y3 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z3 + gene_gpu1.region_pixel * l * unitvector_z;
    							plot(x , y , z , gene_gpu , i , j , k , f_gpu , g_gpu , MU_gpu , &mu_l , &pin_num);
    						}
    		            }
                    }
	            }
	        }
	    }
    }
}

__global__ void rate(float *g_gpu , float *Z_gpu , float *R_gpu , general *gene_gpu){
	general gene_gpu1 = gene_gpu[0];
	int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if(g_gpu[j + i * gene_gpu1.detecter_width + k * gene_gpu1.detecter_width * gene_gpu1.detecter_hight] > 0){
        R_gpu[j + i * gene_gpu1.detecter_width + k * gene_gpu1.detecter_width * gene_gpu1.detecter_hight] = Z_gpu[j + i * gene_gpu1.detecter_width + k * gene_gpu1.detecter_width * gene_gpu1.detecter_hight] / g_gpu[j + i * gene_gpu1.detecter_width + k * gene_gpu1.detecter_width * gene_gpu1.detecter_hight];
    }
    else{
        R_gpu[j + i * gene_gpu1.detecter_width + k * gene_gpu1.detecter_width * gene_gpu1.detecter_hight] = 0;
    }
}

__global__ void mlemreprojection(float *R_gpu , float *h_gpu , general *gene_gpu , colimater *coli_data_gpu){
	general gene_gpu1 = gene_gpu[0];
	int pin_num = 0;
	for(int k = 0; k < gene_gpu1.detecter_num; k++){
		int j = blockIdx.x * blockDim.x + threadIdx.x;//横
	    int i = blockIdx.y * blockDim.y + threadIdx.y;//縦
	    /*検出器の座標(s,t,u座標)*/
		float s = gene_gpu1.radius_rotate + gene_gpu1.distance;
		float t = -(((gene_gpu1.detecter_width - 1.0f) / 2.0f) * gene_gpu1.detecter_pixel) + j * gene_gpu1.detecter_pixel;
		float u = (((gene_gpu1.detecter_hight - 1.0f) / 2.0f) * gene_gpu1.detecter_pixel) - i * gene_gpu1.detecter_pixel;
		float detecter_x = s * cos(coli_data_gpu[k].theta) - t * sin(coli_data_gpu[k].theta);
		float detecter_y = s * sin(coli_data_gpu[k].theta) + t * cos(coli_data_gpu[k].theta);
		float detecter_z = u;
		/*どのコリメータの条件に適しているのかの判定*/
		if(t < coli_data_gpu[k].t1){
            for(pin_num = 0; pin_num <= 2; pin_num++){
                if(pin_num == 0 || pin_num == 2){
                    coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                    coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                    if(pin_num == 0){
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                    }
                    else{
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                    }
                }
                else{
                    coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                    coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                    coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                }
    			float vector_x1 = coli_data_gpu[k].x1 - detecter_x;
    			float vector_y1 = coli_data_gpu[k].y1 - detecter_y;
    			float vector_z1 = coli_data_gpu[k].z1 - detecter_z;
    			float numerator1 = (-vector_x1) * coli_data_gpu[k].housen_x + (-vector_y1) * coli_data_gpu[k].housen_y + (-vector_z1) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
                float denominator1 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x1 , 2.0f) + pow(-vector_y1 , 2.0f) + pow(-vector_z1 , 2.0f))));
                double in_theta1 = acos(numerator1 / denominator1);
                if((in_theta1 * 180.0f / M_PI) <= gene_gpu1.open_degree1){
                    for(int l = 0; l < 1000; l++){
    					/*単位ベクトル*/
    					float vector_size = sqrt(pow(vector_x1 , 2.0f) + pow(vector_y1 , 2.0f) + pow(vector_z1 , 2.0f));
    					float unitvector_x = vector_x1 / vector_size;
    					float unitvector_y = vector_y1 / vector_size;
    					float unitvector_z = vector_z1 / vector_size;
    					float x = coli_data_gpu[k].x1 + gene_gpu1.region_pixel * l * unitvector_x;
    					float y = coli_data_gpu[k].y1 + gene_gpu1.region_pixel * l * unitvector_y;
    					float z = coli_data_gpu[k].z1 + gene_gpu1.region_pixel * l * unitvector_z;
    					replot(x , y , z , gene_gpu , i , j , k , R_gpu , h_gpu , &pin_num);
    				}
                }
            }
        }
		else{
			if(t < coli_data_gpu[k].t2){
                for(pin_num = 0; pin_num <= 6; pin_num++){
                    if(pin_num == 0 || pin_num == 2){
                        coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        if(pin_num == 0){
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                        else{
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                    }
                    else if(pin_num == 1){
                        coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                    }
                    else if(pin_num == 4 || pin_num == 6){
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        if(pin_num == 4){
                            coli_data_gpu[k].z1 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                        else{
                            coli_data_gpu[k].z1 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                        }
                    }
                    else if(pin_num == 5){
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t1 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t1 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                    }
                    else{
                        coli_data_gpu[k].x1 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t1 * sin(coli_data_gpu[k].theta);
                        coli_data_gpu[k].y1 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t1 * cos(coli_data_gpu[k].theta);
                        coli_data_gpu[k].z1 = coli_data_gpu[k].u;
                    }
    				float vector_x1 = coli_data_gpu[k].x1 - detecter_x;
    				float vector_y1 = coli_data_gpu[k].y1 - detecter_y;
    				float vector_z1 = coli_data_gpu[k].z1 - detecter_z;
    				float numerator1 = (-vector_x1) * coli_data_gpu[k].housen_x + (-vector_y1) * coli_data_gpu[k].housen_y + (-vector_z1) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    	            float denominator1 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x1 , 2.0f) + pow(-vector_y1 , 2.0f) + pow(-vector_z1 , 2.0f))));
    	            float in_theta1 = acos(numerator1 / denominator1);
    	            float vector_x2 = coli_data_gpu[k].x2 - detecter_x;
    				float vector_y2 = coli_data_gpu[k].y2 - detecter_y;
    				float vector_z2 = coli_data_gpu[k].z2 - detecter_z;
    				float numerator2 = (-vector_x2) * coli_data_gpu[k].housen_x + (-vector_y2) * coli_data_gpu[k].housen_y + (-vector_z2) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    	            float denominator2 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x2 , 2.0f) + pow(-vector_y2 , 2.0f) + pow(-vector_z2 , 2.0f))));
    	            float in_theta2 = acos(numerator2 / denominator2);
    	            if((in_theta1 * 180.0f / M_PI) <= gene_gpu1.open_degree2){
    	                for(int l = 0; l < 1000; l++){
    						/*単位ベクトル*/
    						float vector_size = sqrt(pow(vector_x1 , 2.0f) + pow(vector_y1 , 2.0f) + pow(vector_z1 , 2.0f));
    						float unitvector_x = vector_x1 / vector_size;
    						float unitvector_y = vector_y1 / vector_size;
    						float unitvector_z = vector_z1 / vector_size;
    						float x = coli_data_gpu[k].x1 + gene_gpu1.region_pixel * l * unitvector_x;
    						float y = coli_data_gpu[k].y1 + gene_gpu1.region_pixel * l * unitvector_y;
    						float z = coli_data_gpu[k].z1 + gene_gpu1.region_pixel * l * unitvector_z;
    						replot(x , y , z , gene_gpu , i , j , k , R_gpu , h_gpu , &pin_num);
    					}
    	            }
    	            else{
    	                if((in_theta2 * 180.0f / M_PI) <= gene_gpu1.open_degree3){
    	                    for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x2 , 2.0f) + pow(vector_y2 , 2.0f) + pow(vector_z2 , 2.0f));
    							float unitvector_x = vector_x2 / vector_size;
    							float unitvector_y = vector_y2 / vector_size;
    							float unitvector_z = vector_z2 / vector_size;
    							float x = coli_data_gpu[k].x2 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y2 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z2 + gene_gpu1.region_pixel * l * unitvector_z;
    							replot(x , y , z , gene_gpu , i , j , k , R_gpu , h_gpu , &pin_num);
    						}
    	                }
    	            }
                }
	        }
	        else{
	            if(t < coli_data_gpu[k].t3){
                    for(pin_num = 0; pin_num <= 6; pin_num++){
                        if(pin_num == 0 || pin_num == 2){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 0){
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 1){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 - (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
                        else if(pin_num == 4 || pin_num == 6){
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 4){
                                coli_data_gpu[k].z2 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z2 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 5){
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t2 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t2 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                        }
                        else{
                            coli_data_gpu[k].x2 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t2 * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y2 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t2 * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z2 = coli_data_gpu[k].u;
                        }
    	                float vector_x2 = coli_data_gpu[k].x2 - detecter_x;
    					float vector_y2 = coli_data_gpu[k].y2 - detecter_y;
    					float vector_z2 = coli_data_gpu[k].z2 - detecter_z;
    					float numerator2 = (-vector_x2) * coli_data_gpu[k].housen_x + (-vector_y2) * coli_data_gpu[k].housen_y + (-vector_z2) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator2 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x2 , 2.0f) + pow(-vector_y2 , 2.0f) + pow(-vector_z2 , 2.0f))));
    	                float in_theta2 = acos(numerator2 / denominator2);
    		            float vector_x3 = coli_data_gpu[k].x3 - detecter_x;
    					float vector_y3 = coli_data_gpu[k].y3 - detecter_y;
    					float vector_z3 = coli_data_gpu[k].z3 - detecter_z;
    					float numerator3 = (-vector_x3) * coli_data_gpu[k].housen_x + (-vector_y3) * coli_data_gpu[k].housen_y + (-vector_z3) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator3 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x3 , 2.0f) + pow(-vector_y3 , 2.0f) + pow(-vector_z3 , 2.0f))));
    		            float in_theta3 = acos(numerator3 / denominator3);
    		            if((in_theta2 * 180.0f / M_PI) <= gene_gpu1.open_degree3){
    		                for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x2 , 2.0f) + pow(vector_y2 , 2.0f) + pow(vector_z2 , 2.0f));
    							float unitvector_x = vector_x2 / vector_size;
    							float unitvector_y = vector_y2 / vector_size;
    							float unitvector_z = vector_z2 / vector_size;
    							float x = coli_data_gpu[k].x2 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y2 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z2 + gene_gpu1.region_pixel * l * unitvector_z;
    							replot(x , y , z , gene_gpu , i , j , k , R_gpu , h_gpu , &pin_num);
    						}
    		            }
    		            else{
    		                if((in_theta3 * 180.0f / M_PI) <= gene_gpu1.open_degree2){
    		                    for(int l = 0; l < 1000; l++){
    								/*単位ベクトル*/
    								float vector_size = sqrt(pow(vector_x3 , 2.0f) + pow(vector_y3 , 2.0f) + pow(vector_z3 , 2.0f));
    								float unitvector_x = vector_x3 / vector_size;
    								float unitvector_y = vector_y3 / vector_size;
    								float unitvector_z = vector_z3 / vector_size;
    								float x = coli_data_gpu[k].x3 + gene_gpu1.region_pixel * l * unitvector_x;
    								float y = coli_data_gpu[k].y3 + gene_gpu1.region_pixel * l * unitvector_y;
    								float z = coli_data_gpu[k].z3 + gene_gpu1.region_pixel * l * unitvector_z;
    								replot(x , y , z , gene_gpu , i , j , k , R_gpu , h_gpu , &pin_num);
    							}
    		                }
    		            }
                    }
	            }
	            else{
                    for(pin_num = 3; pin_num <= 6; pin_num++){
                        if(pin_num == 4 || pin_num == 6){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            if(pin_num == 4){
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u + (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                            else{
                                coli_data_gpu[k].z3 = coli_data_gpu[k].u - (sqrt(1.0f/6.0f) * gene_gpu1.colimater_size);
                            }
                        }
                        else if(pin_num == 5){
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - (coli_data_gpu[k].t3 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + (coli_data_gpu[k].t3 + (sqrt(2.0f/3.0f) * gene_gpu1.colimater_size)) * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
                        else{
                            coli_data_gpu[k].x3 = coli_data_gpu[k].s * cos(coli_data_gpu[k].theta) - coli_data_gpu[k].t3 * sin(coli_data_gpu[k].theta);
                            coli_data_gpu[k].y3 = coli_data_gpu[k].s * sin(coli_data_gpu[k].theta) + coli_data_gpu[k].t3 * cos(coli_data_gpu[k].theta);
                            coli_data_gpu[k].z3 = coli_data_gpu[k].u;
                        }
    	                float vector_x3 = coli_data_gpu[k].x3 - detecter_x;
    					float vector_y3 = coli_data_gpu[k].y3 - detecter_y;
    					float vector_z3 = coli_data_gpu[k].z3 - detecter_z;
    					float numerator3 = (-vector_x3) * coli_data_gpu[k].housen_x + (-vector_y3) * coli_data_gpu[k].housen_y + (-vector_z3) * coli_data_gpu[k].housen_z;//コリメータに入る時の角度を計算するときの分子
    		            float denominator3 = (sqrt(pow(coli_data_gpu[k].housen_x , 2.0f) + pow(coli_data_gpu[k].housen_y , 2.0f) + pow(coli_data_gpu[k].housen_z , 2.0f)) * (sqrt(pow(-vector_x3 , 2.0f) + pow(-vector_y3 , 2.0f) + pow(-vector_z3 , 2.0f))));
    		            float in_theta3 = acos(numerator3 / denominator3);
    		            if((in_theta3 * 180.0f / M_PI) <= gene_gpu1.open_degree1){
    		                for(int l = 0; l < 1000; l++){
    							/*単位ベクトル*/
    							float vector_size = sqrt(pow(vector_x3 , 2.0f) + pow(vector_y3 , 2.0f) + pow(vector_z3 , 2.0f));
    							float unitvector_x = vector_x3 / vector_size;
    							float unitvector_y = vector_y3 / vector_size;
    							float unitvector_z = vector_z3 / vector_size;
    							float x = coli_data_gpu[k].x3 + gene_gpu1.region_pixel * l * unitvector_x;
    							float y = coli_data_gpu[k].y3 + gene_gpu1.region_pixel * l * unitvector_y;
    							float z = coli_data_gpu[k].z3 + gene_gpu1.region_pixel * l * unitvector_z;
    							replot(x , y , z , gene_gpu , i , j , k , R_gpu , h_gpu , &pin_num);
    						}
    		            }
                    }
	            }
	        }
	    }
	}
}

__global__ void divide(float *h_gpu , float *f_gpu , float *C_gpu , general *gene_gpu){
	general gene_gpu1 = gene_gpu[0];
	int j = blockIdx.x * blockDim.x + threadIdx.x;
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    h_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k] = h_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k] / C_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k];
    f_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k] = h_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k] * f_gpu[j + gene_gpu1.image_size * i + gene_gpu1.image_size * gene_gpu1.image_size * k];
}
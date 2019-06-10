#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define H 128
#define W 128


//const float theta = 90.0f * M_PI / 180.0f;
//const char writeFileName[] = "circletouei.raw";
//void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);
void rotateImage(float* f,float* g);
//void copyFloat2Uchar(float* src, unsigned char* dst, const size_t num);

int main()
{
	float f[H*W];
	float g[W*360];
	FILE				*fp1 , *fp2;

    ///////////////////////////////////////////////
    /* 円画像の作成 */
	 for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            //画像の画素(i, j)の座標(x, y)
            float x =  - ( W - 1.0f ) / 2.0f + j;
            float y =    ( H - 1.0f ) / 2.0f - i;
            if(x*x+y*y<50*50){
            	f[H*i+j]=100;
            }
            else{
            	f[H*i+j]=0;
            }

        }
    }
    ///////////////////////////////////////////////

    rotateImage(f,g);
    fp2=fopen("circletouei3.raw","wb");
    fwrite( g ,sizeof( float ) , W * 360 , fp2 );  //8bit
    fclose( fp2 );
    //unsigned char* image = (unsigned char*)calloc(W * 360, sizeof(unsigned char));
    //copyFloat2Uchar(g, image, W * 360);
    
    //writeRawFile(writeFileName, sizeof(unsigned char), W * 360, image);
/************************************************************************************************/
    
   

/************************************************************************************************/
    
    
   
}



void rotateImage(float* f,float* g)
{
	float array[W];
    float theta;
    for(int k=0;k<360;k++)
    {
    	theta=k * M_PI / 180.0f;
    		for (int i = 0; i < H; i++)
    		{
        		for (int j = 0; j < W; j++)
	        		{
	            		//画像の画素(i, j)の座標(x, y)
	            		float x =  - ( W - 1.0f ) / 2.0f + j;
	           			float y =    ( H - 1.0f ) / 2.0f - i;
	            
	            		// 座標(x, y)を原点中心に-theta回転させた座標(X, Y)
	            		float s =   x * cosf(theta) + y * sinf(theta);
                        float t = - x * sinf(theta) + y * cosf(theta);
	            		int idx=floorf(s+W/2);
	            		float dx=s-idx;
	            		if(0<=dx){
	            			g[W*k+idx]+=(1-dx)*f[W*i+j];
	            			g[W*k+idx+1]+=dx*f[W*i+j];
	            		}
	            		if(0>dx){
	            			g[W*k+idx]+=fabs((1-dx))*f[W*i+j];
	            			g[W*k+idx+1]+=fabs(dx)*f[W*i+j];
	            		}
	            		//array[idx] += f[W*i+j];
	            		g[W*k+idx]+=f[W*i+j];
	        		}
	    		
			}

	}
}

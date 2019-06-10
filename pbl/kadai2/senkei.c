#include <stdio.h>
#include <math.h>
void rotateCoordinate(float *f,float *g, float theta);

int main()
{
    float theta = -45.0 / 180.0 * M_PI;
    int             i, j;
    float   f[256*256];       /** Œ´‰æ‘œ—p‚Ì”z—ñ **/
    float   g[256*256];       /** “ñ’l‰»Œã‚Ì‰æ‘œ‚Ì”z—ñ **/
    unsigned char   t[256*256]; 
   
    FILE            *fp1 ,*fp2;

   
	fp1 = fopen("lenna_uchar_256-256.raw" , "rb");
	fread( t , 1 , 256 * 256 , fp1);
	for(i=0;i<256;i++){
	    for(j=0;j<256;j++){
	    	f[256*i+j]=(float)t[256*i+j];
	    }
	}
	fclose( fp1 );

    /****************************************************************/


    /**««ƒvƒƒOƒ‰ƒ€ì¬ŠJŽn««************************************/
	
	rotateCoordinate(f,g,theta);
    /**ªªƒvƒƒOƒ‰ƒ€ì¬I—¹ªª************************************/


    /**ƒgƒŠƒ~ƒ“ƒO‰æ‘œƒtƒ@ƒCƒ‹‘‚«ž‚Ý**************************************/
    fp2 = fopen( "senkei.raw" , "wb" );
    fwrite( g , sizeof(float) , 256 * 256  , fp2 );
    fclose( fp2 );
    /****************************************************************/

    return 0;
}

void rotateCoordinate(float *f,float *g, float theta)
{
    int i,j;
    float X,Y,x,y,x0,y0,a,b;
    for(i=0;i<256;i++){
	    for(j=0;j<256;j++){
	    	x=-(256-1.0f)/2.0f+j;
	    	y=(256-1.0f)/2.0f-i;
	        X=x*cosf(-theta)-y*sinf(-theta);
	        Y=x*sinf(-theta)+y*cosf(-theta);
	        if(-(256.0f-1.0f)/2.0f < X && (256.0f-1.0f)/2.0f > X  && -(256.0f-1.0f)/2.0f < Y && (256.0f-1.0f)/2.0f > Y){
		        x0=floor(X-0.5f)+0.5f;
		        y0=ceil(y-0.5f)+0.5f;
		        a=(256-1)/2-y0;
		        b=x0+(256-1)/2;
		        int index1 = 256*a+b;
		        int index2 = 256*a+b+1;
		        int index3 = 256*(a+1)*b;
		        int index4 = 256*(a+1)*b+1;
		        float dx = fabs(X-x0);
		        float dy = fabs(Y-y0);
		        float s1 = dx*dy;
		        float s2 = (1-dx)*dy;
		        float s3 = dx*(1-dy);
		        float s4 = (1-dx)*(1-dy);
	        	//g[256*i+j]=f[256*a+b]*(1-(X-x0))*(1-(Y-y0))+f[256*a+b+1]*x0*(1-(Y-y0))+f[256*(a+1)+b+1]*x0*y0+f[256*(a+1)+b]*(1-(X-x0))*y0;
	        	//g[256*i+j]=f[256*a+b]*(1-(X-x0))*(1-(y0-Y))+f[256*a+b+1]*(X-x0)*(1-(y0-Y))+f[256*(a+1)+b+1]*(X-x0)*(y0-Y)+f[256*(a+1)+b]*(1-(X-x0))*(y0-Y);
	        	g[256*i+j]=f[index1]*s4 + f[index2]*s3 + f[index3]*s2 + f[index4]*s1;

	        }
	    }
	}
    
}

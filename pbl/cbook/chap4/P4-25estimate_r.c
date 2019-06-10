/* P4-25estimate_r.c */

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#define N 128
#define PI 3.14159265358979

void read_data(char [],float [][N],int); //ファイルからデータを読む
void write_data(char [],float [][N],int); //ファイルにデータを書く
void swap_2D(float[][N],float[][N]); //2次元画像の入れ替え
void dft_2D(int,float[][N],float[][N]); //2次元フーリエ変換
void dft(int,float[],float[]); //1次元フーリエ変換
void phase(float[][N],float[][N],float); //位相画像
void corr(float[][N],float[][N],float[][N],float[][N]); //位相相関
void rotation(float[][N],float[][N],float); //線形補間による画像の回転
float max_data(float[][N]); //最大値

main()
{
	FILE   *fp1;
	char   file_in1[50], file_in2[50];
	float  imr1[N][N]={0.0}, imi1[N][N]={0.0}, imr2[N][N]={0.0};
	float  imr3[N][N]={0.0}, imi3[N][N]={0.0};//配列の初期化
	float  th = (float)0.01;//位相画像の閾値
	int    i, j, ir, k;
	float  deg, r_angle=(float)2.0;//参照画像の回転角度ステップ

	//プログラム題名の表示
	printf("プログラム名：位相相関\n");
	printf("\t2つの画像の位相画像を利用して相関を計算する\n\n");

	printf("観測（平行回転後）画像のファイル名:");
	scanf("%s",file_in1);
	printf("参照（平行回転前）画像のファイル名:");
	scanf("%s",file_in2);

	if((fp1=fopen("corr.txt","w"))==NULL) {
		fprintf(stderr,"fle open error [%s]\n","corr.txt");
		exit(1);
	}

	fprintf(fp1,"k\tcorr\n");

	read_data(file_in1,imr1,N*N);
	read_data(file_in2,imr2,N*N);

	//観測画像の位相画像
	swap_2D(imr1,imi1);
	ir = 1;
	dft_2D(ir,imr1,imi1);

	//画像の位相成分（F(u,v)/|F(u,v)| > th；高周波成分の強調）
	phase(imr1,imi1,th);
	write_data("phase_org.img",imr1,N*N);

	//参照画像の回転と位相画像
	for (k = 0; k < 20; k++) {
		deg = k*r_angle;	
		rotation(imr2,imr3,deg);

		for(i = 0; i < N; i++) {
			for(j = 0; j < N; j++) {
				imi3[i][j]=0.0;
			}
		}

		swap_2D(imr3,imi3);

		ir = 1;

		dft_2D(ir,imr3,imi3);
	
		phase(imr3,imi3,th);

		//観測画像と参照画像の位相相関
		corr(imr1,imi1,imr3,imi3);

		//フーリエ逆変換
		ir = -1;
		dft_2D(ir,imr3,imi3);

		swap_2D(imr3,imi3);

		max_data(imr3);
		printf("k = %d\tmax = %f\n", k,max_data(imr3));
		fprintf(fp1,"%d\t%f\n",k,max_data(imr3));
	}
}

//位相画像
void phase(float a[][N],float b[][N], float r)
{
	float  amp;
	int    i, j;

	for(i = 0; i<N; i++) {
		for(j = 0; j<N; j++) {
			amp = (float)sqrt((double)(a[i][j]*a[i][j] + b[i][j]*b[i][j]));
			if(amp < r) {
				a[i][j] = 0.0;
				b[i][j] = 0.0;
			}
			else
				a[i][j] /= amp;
				b[i][j] /= amp;
		}
	}
}

//フーリエ位相相関（ x1 + iy1 と x2 + iy2 の共役複素数 x2 - iy2 の掛け算）
void corr(float x1[][N],float y1[][N],float x2[][N],float y2[][N])
{
	int i,j;
	float c[N][N]={0.0},d[N][N]={0.0};

	for(i = 0; i<N; i++) {
		for(j = 0; j<N; j++) {
			c[i][j] = x1[i][j]*x2[i][j] + y1[i][j]*y2[i][j];
			d[i][j] = -x1[i][j]*y2[i][j] + y1[i][j]*x2[i][j];
		}
	}

	for(i = 0; i<N; i++) {
		for(j = 0; j<N; j++) {
			x2[i][j] = c[i][j];
			y2[i][j] = d[i][j];
		}
	}
}

//線形補間による画像の回転（反時計回り）
void rotation(float image_in[][N],float image_out[][N], float iv)
{
	int		i,j,x,y,jx,iy;
	float	x0,y0,x1,y1,s,t;
	double	r;
	float	co, si;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			image_out[i][j] = 0.0;
		}
	}

	r=iv*2*PI/360.0;

	co = (float)cos(r);
	si = (float)sin(r);
	for (i = 0; i < N; i++) {	//回転後の配列
		y = N/2 - i;	//回転後のy座標
		for (j = 0; j < N; j++) {
			x = j - N/2;

			s = (x*co + y*si);	//回転前のx座標
			t = (-x*si + y*co);	//回転前のy座標
			
			jx = (int)(s + N/2);	//回転前の配列
			x0 = s + N/2 - jx;
			x1 = 1 - x0;
			iy = (int)(N/2 - t);
			y0 = N/2 - t - iy;
			y1 = 1 - y0;

			if(jx<0||jx>N-2||iy<0||iy>N-2)
			continue;
			
			//線形補間
			image_out[i][j] = y1*x1*image_in[iy][jx] + y0*x1*image_in[iy+1][jx] +
				y1*x0*image_in[iy][jx+1] + y0*x0*image_in[iy+1][jx+1];
		}
	}
}

//相関の最大値を算出
float max_data(float img[][N])
{
	int i,j;
	float max = 0.0;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if(img[i][j] > max) max = img[i][j];
		}
	}
	return max;
}

//2次元画像の入れ替え
void swap_2D(float imr[][N],float imi[][N])
{
//画像の（第1，4象限）と（第2，3象限）の入れ替え

	int i,j,i2,j2;
	float temp;

	for(i=0;i<N;i++){
		for(j=0;j<N/2;j++){
			j2=j+N/2;// 画像中心からj離れた位置を算出

			temp     =imr[i][j];// 実数部についての計算	
			imr[i][j] =imr[i][j2];
			imr[i][j2]=temp;

			temp     =imi[i][j];// 虚数部についての計算	
			imi[i][j] =imi[i][j2];
			imi[i][j2]=temp;

		}
	}

	//画像の（第1，2象限）と（第3，4象限）の入れ替え
	for(j=0;j<N;j++){
		for(i=0;i<N/2;i++){
			i2=i+N/2;// 画像中心からj離れた位置を算出

			temp     =imr[i][j];// 実数部についての計算
			imr[i][j] =imr[i2][j];
			imr[i2][j]=temp;

			temp     =imi[i][j];//虚数部についての計算
			imi[i][j] =imi[i2][j];
			imi[i2][j]=temp;
		}
	}
}

//2次元フーリエ変換
void dft_2D(int ir,float imr[][N],float imi[][N])
{

	float a1[N]={0.0},b1[N]={0.0};
	int i,j;

	for(i = 0; i<N; i++) {
		for(j = 0; j<N; j++) {
			a1[j] = imr[i][j];
			b1[j] = imi[i][j];	
		}		
			dft(ir,a1,b1);//関数に配列の先頭アドレスを渡す

			for(j = 0; j<N; j++) {
			imr[i][j] = a1[j];
			imi[i][j] = b1[j];
			}
	}

	for(j = 0; j<N; j++) {
		for(i = 0; i<N; i++) {
			a1[i] = imr[i][j];
			b1[i] = imi[i][j];
		}

		    dft(ir,a1,b1);
			for(i = 0; i<N; i++) {
				imr[i][j] = a1[i];
				imi[i][j] = b1[i];
			}
	}
}

//1次元フーリエ変換
void dft(int ir,float a[],float b[])
{
	float fr[N]={0.0},fi[N]={0.0};
	int u,x,n;
	
	if(ir == 1) n = 1;
	else
		n = N;

	for(u = 0; u<N; u++) {
		for(x = 0; x<N; x++) {
			fr[u] += (float)(a[x]*cos(2*PI*u*x/N) + ir* b[x]*sin(2*PI*u*x/N));
			fi[u] += (float)(-ir*a[x]*sin(2*PI*u*x/N) + b[x]*cos(2*PI*u*x/N));
		}
	}

	for(x = 0; x<N; x++) {
		a[x] = fr[x]/n;
		b[x] = fi[x]/n;
	}	
}

void read_data(char fi[], float a[][N], int size)
{
	FILE   *fp;

	// open file and read data
	if((fp=fopen(fi,"rb"))==NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fread(a, sizeof(float), size, fp);
	fclose(fp);
}

void write_data(char fi[], float a[][N], int size)
{
	FILE   *fp;

	// open file and write data 
	if((fp=fopen(fi,"wb"))==NULL) {
		fprintf( stderr," Error : file open [%s].\n", fi);
		exit(1);
	}
	fwrite(a, sizeof(float), size, fp);
	fclose(fp);
}

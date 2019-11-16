/* P4-25estimate_r.c */

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#define N 128
#define PI 3.14159265358979

void read_data(char [],float [][N],int); //�t�@�C������f�[�^��ǂ�
void write_data(char [],float [][N],int); //�t�@�C���Ƀf�[�^������
void swap_2D(float[][N],float[][N]); //2�����摜�̓���ւ�
void dft_2D(int,float[][N],float[][N]); //2�����t�[���G�ϊ�
void dft(int,float[],float[]); //1�����t�[���G�ϊ�
void phase(float[][N],float[][N],float); //�ʑ��摜
void corr(float[][N],float[][N],float[][N],float[][N]); //�ʑ�����
void rotation(float[][N],float[][N],float); //���`��Ԃɂ��摜�̉�]
float max_data(float[][N]); //�ő�l

main()
{
	FILE   *fp1;
	char   file_in1[50], file_in2[50];
	float  imr1[N][N]={0.0}, imi1[N][N]={0.0}, imr2[N][N]={0.0};
	float  imr3[N][N]={0.0}, imi3[N][N]={0.0};//�z��̏�����
	float  th = (float)0.01;//�ʑ��摜��臒l
	int    i, j, ir, k;
	float  deg, r_angle=(float)2.0;//�Q�Ɖ摜�̉�]�p�x�X�e�b�v

	//�v���O�����薼�̕\��
	printf("�v���O�������F�ʑ�����\n");
	printf("\t2�̉摜�̈ʑ��摜�𗘗p���đ��ւ��v�Z����\n\n");

	printf("�ϑ��i���s��]��j�摜�̃t�@�C����:");
	scanf("%s",file_in1);
	printf("�Q�Ɓi���s��]�O�j�摜�̃t�@�C����:");
	scanf("%s",file_in2);

	if((fp1=fopen("corr.txt","w"))==NULL) {
		fprintf(stderr,"fle open error [%s]\n","corr.txt");
		exit(1);
	}

	fprintf(fp1,"k\tcorr\n");

	read_data(file_in1,imr1,N*N);
	read_data(file_in2,imr2,N*N);

	//�ϑ��摜�̈ʑ��摜
	swap_2D(imr1,imi1);
	ir = 1;
	dft_2D(ir,imr1,imi1);

	//�摜�̈ʑ������iF(u,v)/|F(u,v)| > th�G�����g�����̋����j
	phase(imr1,imi1,th);
	write_data("phase_org.img",imr1,N*N);

	//�Q�Ɖ摜�̉�]�ƈʑ��摜
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

		//�ϑ��摜�ƎQ�Ɖ摜�̈ʑ�����
		corr(imr1,imi1,imr3,imi3);

		//�t�[���G�t�ϊ�
		ir = -1;
		dft_2D(ir,imr3,imi3);

		swap_2D(imr3,imi3);

		max_data(imr3);
		printf("k = %d\tmax = %f\n", k,max_data(imr3));
		fprintf(fp1,"%d\t%f\n",k,max_data(imr3));
	}
}

//�ʑ��摜
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

//�t�[���G�ʑ����ցi x1 + iy1 �� x2 + iy2 �̋��𕡑f�� x2 - iy2 �̊|���Z�j
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

//���`��Ԃɂ��摜�̉�]�i�����v���j
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
	for (i = 0; i < N; i++) {	//��]��̔z��
		y = N/2 - i;	//��]���y���W
		for (j = 0; j < N; j++) {
			x = j - N/2;

			s = (x*co + y*si);	//��]�O��x���W
			t = (-x*si + y*co);	//��]�O��y���W
			
			jx = (int)(s + N/2);	//��]�O�̔z��
			x0 = s + N/2 - jx;
			x1 = 1 - x0;
			iy = (int)(N/2 - t);
			y0 = N/2 - t - iy;
			y1 = 1 - y0;

			if(jx<0||jx>N-2||iy<0||iy>N-2)
			continue;
			
			//���`���
			image_out[i][j] = y1*x1*image_in[iy][jx] + y0*x1*image_in[iy+1][jx] +
				y1*x0*image_in[iy][jx+1] + y0*x0*image_in[iy+1][jx+1];
		}
	}
}

//���ւ̍ő�l���Z�o
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

//2�����摜�̓���ւ�
void swap_2D(float imr[][N],float imi[][N])
{
//�摜�́i��1�C4�ی��j�Ɓi��2�C3�ی��j�̓���ւ�

	int i,j,i2,j2;
	float temp;

	for(i=0;i<N;i++){
		for(j=0;j<N/2;j++){
			j2=j+N/2;// �摜���S����j���ꂽ�ʒu���Z�o

			temp     =imr[i][j];// �������ɂ��Ă̌v�Z	
			imr[i][j] =imr[i][j2];
			imr[i][j2]=temp;

			temp     =imi[i][j];// �������ɂ��Ă̌v�Z	
			imi[i][j] =imi[i][j2];
			imi[i][j2]=temp;

		}
	}

	//�摜�́i��1�C2�ی��j�Ɓi��3�C4�ی��j�̓���ւ�
	for(j=0;j<N;j++){
		for(i=0;i<N/2;i++){
			i2=i+N/2;// �摜���S����j���ꂽ�ʒu���Z�o

			temp     =imr[i][j];// �������ɂ��Ă̌v�Z
			imr[i][j] =imr[i2][j];
			imr[i2][j]=temp;

			temp     =imi[i][j];//�������ɂ��Ă̌v�Z
			imi[i][j] =imi[i2][j];
			imi[i2][j]=temp;
		}
	}
}

//2�����t�[���G�ϊ�
void dft_2D(int ir,float imr[][N],float imi[][N])
{

	float a1[N]={0.0},b1[N]={0.0};
	int i,j;

	for(i = 0; i<N; i++) {
		for(j = 0; j<N; j++) {
			a1[j] = imr[i][j];
			b1[j] = imi[i][j];	
		}		
			dft(ir,a1,b1);//�֐��ɔz��̐擪�A�h���X��n��

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

//1�����t�[���G�ϊ�
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

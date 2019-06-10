/*符号化したものを10進数にしてないversion*/


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "hahuman.h"
using namespace std;
#define N 8
#define W 256
#define H 256
#define Q 1//量子化のテスト用（通常は1）

const char readfilename[] = "lenna_float_256-256.raw";
const char DCT_image[] 	 = "DCT_float_256-256.raw";
const char result[] 	 = "jpeg_float_256-256.raw";
const char test[] 	 = "code_2_int.raw";

void DCT_quantization(float* f, float* g);
float C(int arg);
void zigzag_scan(float* f, float* g);
void quantization_BrightnessSignal(float* f, float* g);
void quantization_ColorDifferenceSignal(float* f, float* g);
void entropyEncoder_DC(float* f, vector<int> &code, int i, int j, int &dci_1);
void entropyEncoder_AC(float* f, vector<int> &code, int i, int j, int &r);
void encoder(float* f, vector<int> &code);
vector<int> Conversion_10to2(int diff);
void entropyDecoder_DC(vector<int> &code, int &count, float* f, int i, int j, float &dci_1);
float Conversion_2to10(vector<int> code, int &count, int bit_value);
void dequantization_BrightnessSignal(float* f, float* g);
void IDCT(float* f, float* g);
void dezigzag_scan(float* f, float* h);
void entropyDecoder_AC3(vector<int> &code, int &count, float* f, int i, int j, int &r, int &end, float &keep, int &flag);
void decoder(vector<int> &code, float* f);
void do_dezigzag(float* f, float* g);


template <typename T> void 
readRawFile(const char fname[], const size_t num, T* image);

template <typename T> void 
writeRawFile(const char fname[], const size_t num, T* image);


int main(void)
{
	float* f = (float*)calloc(H * W, sizeof(float));
	float* g = (float*)calloc(H * W, sizeof(float));
	float* h1 = (float*)calloc(H * W, sizeof(float));
	float* h2 = (float*)calloc(H * W, sizeof(float));
	float* h3 = (float*)calloc(H * W, sizeof(float));
	float* h4 = (float*)calloc(H * W, sizeof(float));
	vector<int> code;

	readRawFile(readfilename, H * W, f);
	
	DCT_quantization(f, g);

	///////////////////////////////////

	{
		printf("\n");
		int test[8];
		for(int j = N * 0; j < N * 1; j++)
		{
			for(int i = N * 12, k = 0;i < N * 13; i++, k++)
			{
				test[k] = round(g[j * W + i]);
				printf("%d ", test[k]);
			}
			printf("\n");
		}

		printf("\n----this output is in the main function----\n\n");

	}
	///////////////////////////////////

	
	encoder(g, code);






	/////////////////////////////////////////
	/* 符号を2進数のままである場合と10進数に変換した場合の大きさの比較 */
	/*int testarray[62889] = {};
	int count = 0;
	for(int i = 0; i < 62889; i++)
	{
		count++;
		 testarray[i] = code[i];
	}

	cout << "count = " << count << "\n";
	writeRawFile(test, 62889, testarray);
	printf("\n----this output is in the main function----\n\n");*/
	/////////////////////////////////////////







	
	/////////////////////////////////////
	/* 符号語の数と圧縮率*/
	/*int test = 0;
	for(int i = 0; i < code.size(); i++)
	{
		test++;
	}
	printf("code = %d\n", test);
	cout << "Compressibility : " << (float)(test) / (8 * W * H) << "\n";
	printf("---- this printf is in the main function ----\n");*/
	/////////////////////////////////////

	decoder(code, h1);

	///////////////////////////////////
	
	/*{
		printf("\n");
		int test[8];
		for(int j = N * 0; j < N * 1; j++)
		{
			for(int i = N * 2, k = 0;i < N * 3; i++, k++)
			{
				test[k] = round(h1[j * W + i]);
				printf("%d ", test[k]);
			}
			printf("\n");
		}

		printf("\n----this output is in the main function----\n");

	}*/
	///////////////////////////////////


	dezigzag_scan(h1, h2);
	IDCT(h2, h3);
	
	writeRawFile(result, W * H, h3);

	free(f);
	free(g);
	free(h1);
	free(h2);
	free(h3);
	free(h4);
}


void dezigzag_scan(float* f, float* h)
{

	float* g = (float*)calloc(N * N, sizeof(float));
	float* h0 = (float*)calloc(N * N, sizeof(float));
	float* h1 = (float*)calloc(N * N, sizeof(float));

	for(int m = 0; m < H / N; m++)
	{
		for(int n = 0; n < W / N; n++)
		{
			for(int j0 = N * m, j1 = 0; j1 < N; j0++, j1++)
			{
				for(int k0 = N * n, k1 = 0; k1 < N; k0++, k1++)
				{
					g[j1 * N + k1] = f[j0 * W + k0];
				}
			}
			do_dezigzag(g, h0);
			dequantization_BrightnessSignal(h0, h1);

			for(int j = N * m, u = 0; u < N; j++, u++)
			{
				for(int k = N * n, v = 0; v < N; k++, v++)
				{
					h[j * W + k] = h1[u * N + v];					
				}
			}
		}
	}

	
	free(g);
	free(h0);
	free(h1);

}


void IDCT(float* f, float* g)
{
	for(int m = 0; m < H / N; m++)
	{
		for(int n = 0; n < W / N; n++)
		{
			for(int j0 = N * m, j1 = 0; j1 < N; j0++, j1++)
			{
				for(int k0 = N * n, k1 = 0; k1 < N; k0++, k1++)
				{
					float sum = 0.0f;
					for(int u = 0, u1 = N * m; u < N; u++, u1++)
					{
						for(int v = 0, v1 = N * n; v < N; v++, v1++)
						{
							sum += C(u) * C(v) * f[u1 * W + v1] * cos(((2 * j1 + 1) * u * M_PI) / 16) * cos(((2 * k1 + 1) * v * M_PI) / 16);
						}
					}
					g[j0 * W + k0] = sum / 4;
				}
			}

			/*for(int j = N * m, u = 0; u < N; j++, u++)
			{
				for(int k = N * n, v = 0; v < N; k++, v++)
				{
					//DCTの計算結果を画像として出力する場合
					//i[j * W + k] = g0[u * N + v];
					//ジグザグスキャンまでの結果を出力する場合
					i[j * W + k] = g2[u * N + v];					
					g0[u * N + v] = 0.0f;
				}
			}*/
		}
	}


}

void dequantization_BrightnessSignal(float* f, float* g)
{
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			g[i * N + j] = f[i * N + j] * (q_table[i][j] * Q);
		}
	}
}

float Conversion_2to10(vector<int> code, int &count, int bit_value)
{
	int minus = 0;

	if(code[count] == 0)
	{
		minus = 1;
	}

	//補数表現をもとにもどす
	if(minus)
	{
		for(int i = 0; i < bit_value; i++)
		{
			if(code[count + i] == 1)
				code[count + i] = 0;
			else
				code[count + i] = 1;
		}
	}

	float result = 0.0f;
	for(int i = bit_value - 1, j = count; i >= 0; i--, j++)
	{
		result += code[j] * pow(2.0f, i);
	}

	if(minus)
		result *= -1;

	count += bit_value;

	return result;
}

void entropyDecoder_DC(vector<int> &code, int &count, float* f, int i, int j, float &dci_1)
{

	int bit_value;
	int diff;
	int ssss;
	/* dc_code_tableから一致する符号を調べる */
	for(int j = 0; j < 12; j++)
	{
		int count_start = count;
		for(int k = 0; k < 9; k++)
		{
			if(dc_code_table[j][k] == code[count])
			{
				//もしdc_code_table[j][k + 1]が-1ならそのjがbit_valueになる
				if(dc_code_table[j][k + 1] == -1)
				{			
					bit_value = j;
					k = 9;
					j = 12;//完全一致のため調べるのは終了
				}
				count++;
			}
			else
			{
				//もし一致しなかったら次の符号表と比較
				k = 9;
				count = count_start;//違った場合、元の場所までもどす
			}
		}
	}


	/*bit_valueビット2進10進変換*/
	diff = Conversion_2to10(code, count, bit_value);
	f[i * W + j] = diff + dci_1;
	dci_1 = f[i * W + j];
}


void do_dezigzag(float* f, float* g)
{
	g[0 * N + 0] = f[0 * N + 0];
	g[0 * N + 1] = f[0 * N + 1];
	g[1 * N + 0] = f[0 * N + 2];
	g[2 * N + 0] = f[0 * N + 3];
	g[1 * N + 1] = f[0 * N + 4];
	g[0 * N + 2] = f[0 * N + 5];
	g[0 * N + 3] = f[0 * N + 6];
	g[1 * N + 2] = f[0 * N + 7];

	g[2 * N + 1] = f[1 * N + 0];
	g[3 * N + 0] = f[1 * N + 1];
	g[4 * N + 0] = f[1 * N + 2];
	g[3 * N + 1] = f[1 * N + 3];
	g[2 * N + 2] = f[1 * N + 4];
	g[1 * N + 3] = f[1 * N + 5];
	g[0 * N + 4] = f[1 * N + 6];
	g[0 * N + 5] = f[1 * N + 7];

	g[1 * N + 4] = f[2 * N + 0];
	g[2 * N + 3] = f[2 * N + 1];
	g[3 * N + 2] = f[2 * N + 2];
	g[4 * N + 1] = f[2 * N + 3];
	g[5 * N + 0] = f[2 * N + 4];
	g[6 * N + 0] = f[2 * N + 5];
	g[5 * N + 1] = f[2 * N + 6];
	g[4 * N + 2] = f[2 * N + 7];
	
	g[3 * N + 3] = f[3 * N + 0];
	g[2 * N + 4] = f[3 * N + 1];
	g[1 * N + 5] = f[3 * N + 2];
	g[0 * N + 6] = f[3 * N + 3];
	g[0 * N + 7] = f[3 * N + 4];
	g[1 * N + 6] = f[3 * N + 5];
	g[2 * N + 5] = f[3 * N + 6];
	g[3 * N + 4] = f[3 * N + 7];

	g[4 * N + 3] = f[4 * N + 0];
	g[5 * N + 2] = f[4 * N + 1];
	g[6 * N + 1] = f[4 * N + 2];
	g[7 * N + 0] = f[4 * N + 3];
	g[7 * N + 1] = f[4 * N + 4];
	g[6 * N + 2] = f[4 * N + 5];
	g[5 * N + 3] = f[4 * N + 6];
	g[4 * N + 4] = f[4 * N + 7];

	g[3 * N + 5] = f[5 * N + 0];
	g[2 * N + 6] = f[5 * N + 1];
	g[1 * N + 7] = f[5 * N + 2];
	g[2 * N + 7] = f[5 * N + 3];
	g[3 * N + 6] = f[5 * N + 4];
	g[4 * N + 5] = f[5 * N + 5];
	g[5 * N + 4] = f[5 * N + 6];
	g[6 * N + 3] = f[5 * N + 7];

	g[7 * N + 2] = f[6 * N + 0];
	g[7 * N + 3] = f[6 * N + 1];
	g[6 * N + 4] = f[6 * N + 2];
	g[5 * N + 5] = f[6 * N + 3];
	g[4 * N + 6] = f[6 * N + 4];
	g[3 * N + 7] = f[6 * N + 5];
	g[4 * N + 7] = f[6 * N + 6];
	g[5 * N + 6] = f[6 * N + 7];

	g[6 * N + 5] = f[7 * N + 0];
	g[7 * N + 4] = f[7 * N + 1];
	g[7 * N + 5] = f[7 * N + 2];
	g[6 * N + 6] = f[7 * N + 3];
	g[5 * N + 7] = f[7 * N + 4];
	g[6 * N + 7] = f[7 * N + 5];
	g[7 * N + 6] = f[7 * N + 6];
	g[7 * N + 7] = f[7 * N + 7];




}



void DCT_quantization(float* f, float* i)
{
	float* g0 = (float*)calloc(N * N, sizeof(float));
	float* g1 = (float*)calloc(N * N, sizeof(float));
	float* g2 = (float*)calloc(N * N, sizeof(float));

	for(int m = 0; m < H / N; m++)
	{
		for(int n = 0; n < W / N; n++)
		{
			for(int u = 0; u < N; u++)
			{
				for(int v = 0; v < N; v++)
				{
					float sum = 0.0f;
					for(int j0 = N * m, j1 = 0; j1 < N; j0++, j1++)
					{
						for(int k0 = N * n, k1 = 0; k1 < N; k0++, k1++)
						{
							sum += f[j0 * W + k0] * cos(((2 * j1 + 1) * u * M_PI) / 16) * cos(((2 * k1 + 1) * v * M_PI) / 16);
						}
					}
					g0[u * N + v] = (C(u) * C(v) / 4) * sum;
				}
			}
			quantization_BrightnessSignal(g0, g1);
			zigzag_scan(g1, g2);

			for(int j = N * m, u = 0; u < N; j++, u++)
			{
				for(int k = N * n, v = 0; v < N; k++, v++)
				{
					//DCTの計算結果を画像として出力する場合
					//i[j * W + k] = g0[u * N + v];
					//ジグザグスキャンまでの結果を出力する場合
					i[j * W + k] = g2[u * N + v];					
					g0[u * N + v] = 0.0f;
				}
			}
		}
	}
	free(g0);
	free(g1);
	free(g2);

}

float C(int arg)
{
	if(arg == 0)
		return 1 / sqrt(2);
	else
		return 1;
}

void zigzag_scan(float* f, float* g)
{
	g[0] = f[0];
	g[1] = f[0 * N + 1];
	g[2] = f[1 * N + 0];
	{			
		int count = 3;
		int i = 2;
		int j = 0;
		int i_max = 2;
		int j_max = 3;
		int i_up = 1;
		int j_up = 1;
		while(count != 36)
		{
			g[count] = f[i * N + j];
			//printf("count = %d	j = %d	i = %d\n", count, j, i);
			if(i_up)
			{
				if(i == i_max)
				{
					i--;
					i_max += 2;
					i_up = 0;
				}
				else
					i++;
			}
			else
			{
				if(i == 0)
					i_up = 1;
				else
					i--;
			}

			if(j_up)
			{
				if(j_max == j)
				{
					j--;
					j_max += 2;
					j_up = 0;
				}
				else
					j++;
			}
			else
			{
				if(j == 0)
					j_up = 1;
				else
					j--;
			}

			count++;
		}

		int i_min = 1;
		int j_min = 2;
		j = 1;
		i = 7;
		g[36] = f[7 * N + 1];
		//printf("count = %d	j = %d	i = %d\n", count, j, i);
		count++;
		i--;
		j++;
		i_up = 0;

		while(count != N * N)
		{
			i_max = 7;
			j_max = 7;
			g[count] = f[i * N + j];
			//printf("count = %d	j = %d	i = %d\n", count, j, i);
			if(i_up)
			{
				if(i == i_max)
				{
					i_up = 0;
				}
				else
					i++;
			}
			else
			{
				if(i == i_min)
				{
					i_min += 2;
					i++;
					i_up = 1;
				}
				else
					i--;
			}

			if(j_up)
			{
				if(j_max == j)
				{
					j_up = 0;
				}
				else
					j++;
			}
			else
			{
				if(j == j_min)
				{
					j_min += 2;
					j++;
					j_up = 1;
				}
				else
					j--;
			}

			count++;

		}
	}
}

void quantization_BrightnessSignal(float* f, float* g)
{


	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			g[i * N + j] = f[i * N + j] / (q_table[i][j] * Q);//Qはテスト用
	}
}


void quantization_ColorDifferenceSignal(float* f, float* g)
{
	int table[N * N] = {	17, 18, 24, 47, 99, 99, 99, 99,
		  			    	18, 21, 26, 66, 99, 99, 99, 99,
			    		 	24, 26, 56, 66, 99, 99, 99, 99,
			    		 	47, 66, 99, 99, 99, 99, 99, 99,
				         	99, 99, 99, 99, 99, 99, 99, 99,
				    	 	99, 99, 99, 99, 99, 99, 99, 99,
				    		99, 99, 99, 99, 99, 99, 99, 99,
				    		99, 99, 99, 99, 99, 99, 99, 99, };

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			g[i * N + j] = f[i * N + j] / table[i * N + j];
	}

}


void encoder(float* f, vector<int> &code)
{
	int r = 0;
	int dci_1 = 0;

	for(int m = 0; m < H / N; m++)
	{
		for(int n = 0; n < W / N; n++)
		{
			for(int j0 = N * m, j1 = 0; j1 < N; j0++, j1++)
			{
				for(int k0 = N * n, k1 = 0; k1 < N; k0++, k1++)
				{
					if(j1 == 0 && k1 == 0)
					{
						entropyEncoder_DC(f, code, j0, k0, dci_1);
					}
					else
						entropyEncoder_AC(f, code, j0, k0, r);
				}
			}
		}
	}


	/////////////////////////////////////////////
	/* test */
	//printf("code result\n\n");
	//for(int m = 0; m < code.size(); m++)
	/*for(int m = 0; m <50;m++)
	{
		printf("code[%03d] = %d\n", m, code[m]);
	}
	printf("----this output is in the encoder function----\n");*/
	/////////////////////////////////////////////


}


void entropyEncoder_DC(float* f, vector<int> &code, int i, int j, int &dci_1)
{
	int ssss;
	int dci;
	
	dci = round(f[i * W + j]);


	int diff = dci - dci_1;
	//printf("%d\n%f\n", diff, f[i * W + j]);
	//diffが0の時にsegmentation faultになるかもしれない(if文を入れてエラーを回避)
	if(diff != 0)
		ssss = log(abs(diff)) / log(2) + 1;
	else
		ssss = 0;

	/////////////////////////////////////////////
	/* test */
	//printf("dci = %d dci_1 = %d ssss = %d\n", dci, dci_1, ssss);
	/////////////////////////////////////////////

	/*符号語を代入*/
	for(int k = 0; k < dc_length_table[ssss]; k++)
	{
		code.push_back(dc_code_table[ssss][k]);
	}

	if(ssss != 0)
	{
		vector<int> binary;
		binary = Conversion_10to2(diff);
		
		/*逆順で書き込み*/
		for(int k = binary.size() - 1; k >= 0; k--)
		{
			code.push_back(binary[k]);
		}
	}
	dci_1 = dci;
}



void entropyEncoder_AC(float* f, vector<int> &code, int i, int j, int &r)
{
	int ac = 0;

	ac = round(f[i * W + j]);

	/*AC成分が正のとき、AC成分を2進数に変換して書き込み*/
	/*AC成分が負のとき、1の補数表現をとって書き込み*/
	if(ac != 0)
	{
		int ssss = log(abs(ac)) / log(2) + 1;

		for(int k = 0; k < ac_length_table[r * 11 + ssss]; k++)
		{
			code.push_back(ac_code_table[r * 11 + ssss][k]);
		}
		vector<int> binary;
		binary = Conversion_10to2(ac);
		
		/*逆順で書き込み*/
		for(int k = binary.size() - 1; k >= 0; k--)
		{
			code.push_back(binary[k]);
		}	

		r = 0;
	}
	else
	{
		//もしac = 0だった場合ゼロランを1追加(maxで15)
		if(r < 15)
			r++;
	}
	//ブロックの最後であれば1010(EOB)を入れる
	if((i + 1) % 8 == 0 && (j + 1) % 8 == 0)
	{
		for(int k = 0; k < ac_length_table[0]; k++)
		{
			code.push_back(ac_code_table[0][k]);
		}

		r = 0;
	}

}


vector<int> Conversion_10to2(int diff)
{
	vector<int> v;
	int minus = 0;

	if(diff < 0)
		minus = 1;
	
	diff = abs(diff);

	for(int i = diff; i > 0; i /= 2)
		v.push_back(i % 2);

	if(minus)
	{
		/*DIFFの１の補数表現をとる*/
		for(int k = 0; k < v.size(); k++)
		{
			if(v[k] == 1)
				v[k] = 0;
			else
				v[k] = 1;
		}
	}

	return v;

}


void decoder(vector<int> &code, float* f)
{
	int count = 0;
	int r = 0;
	int flag = 0;
	float keep = 0.0f;
	float dci_1 = 0.0f;
	//printf("code.size() = %d\n" , (int)code.size());
	for(int m = 0; m < H / N; m++)
	{
		for(int n = 0; n < W / N; n++)
		{
			int end = 0;
			for(int j0 = N * m, j1 = 0; j1 < N; j0++, j1++)
			{
				for(int k0 = N * n, k1 = 0; k1 < N; k0++, k1++)
				{
					if(j1 == 0 && k1 == 0)
					{
						entropyDecoder_DC(code, count, f, j0, k0, dci_1);
					}
					else if(r > 0)
					{
						f[j0 * W + k0] = 0.0f;
						r--;
						flag = 1;
					}
					else if(r == 0 && flag == 1)
					{
						flag = 0;
						f[j0 * W + k0] = keep;
					}

					else if(r == 0 && flag == 0)
						entropyDecoder_AC3(code, count, f, j0, k0, r, end, keep, flag);
				}
			}
		}
	}
	//printf("count = %d\n", count);
	//printf("----this printf is in the decoder function-----\n");
}


void entropyDecoder_AC3(vector<int> &code, int &count, float* f, int i, int j, int &r, int &end, float &keep, int &flag)
{
	int ssss;

	if(!end)
	{
		/* ac_code_tableから一致する符号を調べる */
		for(int m = 0; m < 176; m++)
		{
			int count_start = count;
			for(int n = 0; n < 16; n++)
			{
				if(ac_code_table[m][n] == code[count])
				{
					if(n == 15)//r = 15 ssss = 1のとき、後ろに-1が存在しないための対策
					{			
						ssss = m % 11;
						r = m / 11;
						m = 176;
						n =  16;//完全一致のため調べるのは終了
						//もしEOB(r = 0, ssss = 0)のときブロックの残りは全部
						if(r == 0 && ssss == 0)
						{
							end = 1;
						}
					}
					else
					{
						//もしac_code_table[m][n + 1]が-1ならそのmがssssになる
						if(ac_code_table[m][n + 1] == -1)
						{			
							ssss = m % 11;
							r = m / 11;
							m = 176;
							n =  16;//完全一致のため調べるのは終了
							//もしEOB(r = 0, ssss = 0)のときブロックの残りは全部
							if(r == 0 && ssss == 0)
							{
								end = 1;
							}
						}
					}
					count++;
				}
				else
				{
					//もし一致しなかったら次の符号表と比較
					n = 16;
					count = count_start;
				}
			}
		}

		if(r == 0 && ssss == 0)
		{
			end = 1;
		}
		else if (r == 0)
		{
			f[i * W + j] = Conversion_2to10(code, count, ssss);
		}
		else
		{
			keep = Conversion_2to10(code, count, ssss);
			f[i * W + j] = 0.0f;
			r--;
			if(r == 0)
				flag = 1;
		}

 	}
 	else if(end)
 	{
 		f[i * W + j] = 0.0f;
 	}

}




template <typename T> void 
readRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"rb");

	if(fp == NULL)
	{
		printf("failed to open %s.\n",fname);
		exit(-1);
	}

	size_t ret = fread(image, sizeof(T), num,fp);

	if(num != ret)
	{
		printf("failed to read %s\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}


template <typename T> void 
writeRawFile(const char fname[], const size_t num, T* image)
{
	FILE* fp = fopen(fname,"wb");

	if(fp == NULL)
	{
		printf("failed to open %s\n", fname);
		exit(-1);
	}

	size_t ret = fwrite(image, sizeof(T), num, fp);

	if(num != ret)
	{
		printf("failed to write %s.\n", fname);
		fclose(fp);
		exit(-1);
	}

	fclose(fp);
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 5

typedef struct
{
	char	name[20];
	int		height;
	float	weight;

}gstudent;

void swap_student(gstudent *x, gstudent *y);
void sort_by_height(gstudent *student);

void main(void)
{
	FILE	*fp1, *fp2;
	gstudent *student;	

	student = (gstudent*)malloc(sizeof(gstudent) * N);
	if(student == NULL)
	{
		printf("failed to allocate memory of student.\n");
		exit(1);
	}
	fp1 = fopen("struct.txt","r");
	if(fp1 == NULL)
	{
		printf("failed to open file.\n");
		exit(1);
	}

	for(int i = 0;i < N; i++)
		fscanf(fp1, "%s %d %f", &(student + i) -> name, &(student + i) -> height, &(student + i) -> weight);
	
	fclose(fp1);

	sort_by_height(student);

	fp2 = fopen("struct_new.txt", "w");
	if(fp2 == NULL)
	{
		printf("failed to open file.\n");
		exit(1);
	}

	for(int i = 0;i < N; i++)
		fprintf(fp2, "%6s %3d %2.1f\n", (student + i) -> name, (student + i) -> height, (student + i) -> weight);

	fclose(fp2);

	free(student);
}

void swap_student(gstudent *x, gstudent *y)
{
	gstudent temp;

	temp = *x;
	*x = *y;
	*y = temp;
}


void sort_by_height(gstudent *student)
{
	for(int i = 0; i < N; i++)
	{
		for(int j = N - 1; j > i; j--)
		{
			if((student + j - 1) -> height > (student + j) -> height)
				swap_student(student + j - 1, student + j);
		}
	}
}
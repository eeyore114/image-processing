#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 2



typedef struct
{
	char	name[20];
	int		height;
	float	weight;
	long	age;
} gstudent;


void function(gstudent *ogawa)
{
	strcpy(ogawa[0].name,"baba");
	ogawa[0].height=3;
	ogawa[0].weight=2.0;
	ogawa[0].age=3;
	
	printf("......\nname=%s\n",ogawa[0].name);
	printf("height=%d\n",ogawa[0].height);
	printf("weight=%f\n",ogawa[0].weight);
	printf("age=%ld\n\n\n",ogawa[0].age);		

	strcpy((ogawa+1)->name,"Yamada");
	(ogawa+1)->height=11;
	(ogawa+1)->weight=11.1;
	(ogawa+1)->age=111;

	printf("->->->->->->->\nname=%s\n",(ogawa+1)->name);
	printf("height=%d\n",(ogawa+1)->height);
	printf("weight=%f\n",(ogawa+1)->weight);
	printf("age=%ld\n",(ogawa+1)->age);	
}


int main()
{		
	gstudent *ogawa;
	ogawa = (gstudent*)malloc(sizeof(gstudent) * N);		
		
	function(ogawa);
		
    return 0;
}


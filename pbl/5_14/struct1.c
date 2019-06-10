#include <stdio.h>
#include <string.h>



struct gstudent
{
	char	name[20];
	int 	height;
	float 	weight;
	long 	age;
};


int main()
{
	struct gstudent a;
	strcpy(a.name,"baba");
	a.height=3;
	a.weight=2.0;
	a.age=3;
		
	printf("name=%s\n",a.name);
	printf("height=%d\n",a.height);
	printf("weight=%f\n",a.weight);
	printf("age=%d\n",a.age);
		
    return 0;
}


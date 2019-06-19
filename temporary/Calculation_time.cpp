#include <iostream>
#include <time.h>
#include <unistd.h>

void start_time(clock_t &start);
void end_time(clock_t start);

int main()
{
	clock_t start;
	start_time(start);
	sleep(5);

	end_time(start);
}

void start_time(clock_t &start)
{
	start = clock();
}
void end_time(clock_t start)
{
	clock_t end = clock();
	float time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
	std::cout << "time = " << time << "[s]\n"
  					 <<	"     = " << time / 60 << "[min]\n"
  					 << "     = " << time / 3600 << "[h]" << std::endl;
}

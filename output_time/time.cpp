#include <stdio.h>
#include <iostream>
#include <time.h>
	
int main(void)
{
    time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */

    strftime(str, 255, "%Y %B %d %A %p %I:%M:%S", date);
    std::string s = str;
    return s;
}
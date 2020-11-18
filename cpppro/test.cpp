#include<stdio.h>
#include<stdlib.h>
#include<string.h>
int main(int argc,char* argv[]){


  char *p;
  p = strtok(argv[1], " ");
  double a = atof(p);
  printf("%lf\n",a);
  p = strtok(NULL, " ");
  printf("%s\n",p);
  float b = atof(p);
  printf("%f\n",b);
  p = strtok(NULL, " ");
  if(p==NULL)
    printf("yes");
  p = strtok(NULL, " ");
  double c = atof(p);
  if(p==NULL)
    printf("%lf",c);
  return 0;
}



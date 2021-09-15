# include<iostream>
# include<fstream>
# include<cmath>

using namespace std;
// to find maximum in an array
double maximum(double A[], int size)
{
	int i;
	double t;
	t=A[1];
	for(i=2;i<size-1;i++)
    {
		if(A[i]>t)
		t=A[i];
	}
	return(t);
}
//  


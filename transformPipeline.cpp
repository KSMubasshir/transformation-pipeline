#include <fstream>
#include <iostream>
#include <stack>
#include <vector>
#include <cmath>
#include <iomanip>
#define pi (2*acos(0.0))
using namespace std;

inline double degreeToRadian(int degree){
    return degree * pi/180;
}

typedef std::vector<std::vector<double> > two_d_vector;
double M[4][4];
double trMatrix[4][4];
double temp[4][4];
double triangle[4][4];
stack<two_d_vector> mystack;
two_d_vector tos;


//==============================Modeling Transformation Functions ==============================

void init(){
	for(int i=0;i<4;i++){
		for (int j = 0; j < 4; j++)
		{
			if(i==j)
				trMatrix[i][j] = 1 ;
			else
				trMatrix[i][j] = 0 ;
		}
	}

	for (int i = 0; i < 4; ++i){	
		std::vector<double> v;
		tos.push_back(v);
		for (int j = 0; j < 4; ++j)
		{	
			if(i==j)
				tos[i].push_back(1); 
			else
				tos[i].push_back(0);
		}
	}
	mystack.push(tos);
}


struct point{
	double x;
	double y;
	double z;
	double w;
};
typedef point Point;

struct glLookAt{
	Point p[3];
	double fov;
	double aspectRatio;
	double near;
	double far;
};
struct glLookAt gLAt;

void matrixMultiplication(){
	for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            temp[i][j] = 0;
            for (int k = 0; k < 4; k++)
            {
                temp[i][j] += tos[i][k] * M[k][j];
            }
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            tos[i][j] = temp[i][j];
        }
    }
}

void handleTriangle(){
	tos = mystack.top();

	for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            temp[i][j] = 0;
            for (int k = 0; k < 4; k++)
            {
                temp[i][j] +=  tos[i][k] * triangle[k][j]  ;
            }
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
        	if(temp[3][i] == 1 )
        		triangle[j][i] = temp[j][i];
        	else //if(temp[3][i] != 0)
        		triangle[j][i] = temp[j][i]/(temp[3][i]*1.0);
        }
    }
}

void handleTranslation(double tx,double ty,double tz){
	M[0][0]=1; M[0][1]=0; M[0][2]=0; M[0][3]=tx;
	M[1][0]=0; M[1][1]=1; M[1][2]=0; M[1][3]=ty;
	M[2][0]=0; M[2][1]=0; M[2][2]=1; M[2][3]=tz;
	M[3][0]=0; M[3][1]=0; M[3][2]=0; M[3][3]=1;
	tos = mystack.top();
	mystack.pop();
	matrixMultiplication();
	mystack.push(tos);
}

void handleRotation(double angle,double ax,double ay,double az){
	double cosTheta = cos(degreeToRadian(angle));
	double sinTheta = sin(degreeToRadian(angle));
	Point a,i,j,k,c1,c2,c3;
	double mod_a = 0;

	mod_a = sqrt( ax*ax + ay*ay + az*az );
	ax /= mod_a ; ay /= mod_a ; az /= mod_a ; //a.normalize

	a.x = ax ; a.y = ay ; a.z = az ; a.w = 0  ;
	i.x = 1 ; i.y = 0 ; i.z = 0 ; i.w = 0  ;
	j.x = 0 ; j.y = 1 ; j.z = 0 ; j.w = 0  ;
	k.x = 0 ; k.y = 1 ; k.z = 0 ; k.w = 0  ;

	c1.x = (cosTheta*i.x) + (1 - cosTheta)*(a.x*a.x) ; 
	c1.y = (1 - cosTheta)*(a.x*a.y) +  sinTheta * (i.x * a.z) ; 
	c1.z = (1 - cosTheta)*(a.x*a.z) -  sinTheta * (i.x * a.y) ; 
	c1.w = 0;

	c2.x = (1 - cosTheta)*(a.y*a.x) - sinTheta * ( j.y * a.z) ; 
	c2.y = (cosTheta*j.y) + (1 - cosTheta)*(a.y*a.y); 
	c2.z = (1 - cosTheta)*(a.y*a.z) +  sinTheta * (a.x * j.y); 
	c2.w = 0 ;

	c3.x = (1 - cosTheta)*(a.z*a.x) + sinTheta * ( a.y * k.z); 
	c3.y = (1 - cosTheta)*(a.z*a.y) - sinTheta * ( a.x * k.z); 
	c3.z = (cosTheta*k.z) + (1 - cosTheta)*(a.z*a.z);; 
	c3.w = 0 ;

	M[0][0]=c1.x;  M[0][1]=c2.x;  M[0][2]=c3.x;  M[0][3]=0;
	M[1][0]=c1.y;  M[1][1]=c2.y;  M[1][2]=c3.y;  M[1][3]=0;
	M[2][0]=c1.z;  M[2][1]=c2.z;  M[2][2]=c3.z;  M[2][3]=0;
	M[3][0]=0;     M[3][1]=0;     M[3][2]=0;     M[3][3]=1;

	tos = mystack.top();
	mystack.pop();
	matrixMultiplication();
	mystack.push(tos);
}

void handleScaling(double sx,double sy,double sz){
	M[0][0]=sx; M[0][1]=0;  M[0][2]=0;  M[0][3]=0;
	M[1][0]=0;  M[1][1]=sy; M[1][2]=0;  M[1][3]=0;
	M[2][0]=0;  M[2][1]=0;  M[2][2]=sz; M[2][3]=0;
	M[3][0]=0;  M[3][1]=0;  M[3][2]=0;  M[3][3]=1;
	tos = mystack.top();
	mystack.pop();
	matrixMultiplication();
	mystack.push(tos);
}

void handlePush(){
	tos = mystack.top();
	mystack.push(tos);
}

void handlePop(){
	mystack.pop();
}


//==============================View Transformation Functions ==============================

double V[4][4];
void calculateViewTransformationMatrix(){
	Point l,r,u;
	double R[4][4];
	double T[4][4];

	l.x = gLAt.p[1].x - gLAt.p[0].x ; 
	l.y = gLAt.p[1].y - gLAt.p[0].y ; 
	l.z = gLAt.p[1].z - gLAt.p[0].z ;

	double mod_l = 0;

	mod_l = sqrt( l.x*l.x + l.y*l.y + l.z*l.z );
	l.x /= (mod_l*1.0) ; l.y /= (mod_l*1.0) ; l.z /= (mod_l*1.0) ; //l.normalize

	//r = l x up
	r.x = l.y * gLAt.p[2].z - gLAt.p[2].y * l.z;
    r.y = gLAt.p[2].x * l.z - l.x * gLAt.p[2].z;
    r.z = l.x * gLAt.p[2].y - gLAt.p[2].x * l.y;

    double mod_r = 0;

	mod_r = sqrt( r.x*r.x + r.y*r.y + r.z*r.z );
	r.x /= (mod_r*1.0) ; r.y /= (mod_r*1.0) ; r.z /= (mod_r*1.0) ; //r.normalize

	//u = r x l
	u.x = r.y * l.z - l.y * r.z;
    u.y = l.x * r.z - r.x * l.z;
    u.z = r.x * l.y - l.x * r.y;


    T[0][0]=1;  T[0][1]=0;  T[0][2]=0;  T[0][3]= - gLAt.p[0].x;
	T[1][0]=0;  T[1][1]=1;  T[1][2]=0;  T[1][3]= - gLAt.p[0].y;
	T[2][0]=0;  T[2][1]=0;  T[2][2]=1;  T[2][3]= - gLAt.p[0].z;
	T[3][0]=0;  T[3][1]=0;  T[3][2]=0;  T[3][3]= 1 ;

    R[0][0]=  r.x;  R[0][1]=  r.y;  R[0][2]=  r.z;  R[0][3]=0;
	R[1][0]=  u.x;  R[1][1]=  u.y;  R[1][2]=  u.z;  R[1][3]=0;
	R[2][0]= -l.x;  R[2][1]= -l.y;  R[2][2]= -l.z;  R[2][3]=0;
	R[3][0]=  0;    R[3][1]=  0;    R[3][2]=  0;    R[3][3]=1;


	for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            V[i][j] = 0;
            for (int k = 0; k < 4; k++)
            {
                V[i][j] += R[i][k] * T[k][j];
            }
        }
    }
}


//==============================Projection Transformation Functions ==============================

double P[4][4];
void calculateProjectionTransformationMatrix(){
	double fovY =  gLAt.fov ;
	double fovX =  fovY  *  gLAt.aspectRatio ;
	double near =  gLAt.near ;
	double far =  gLAt.far ;
	double t = near * tan(degreeToRadian(fovY/2.0)) ;
	double r = near * tan(degreeToRadian(fovX/2.0)) ;

	P[0][0]=near/r ;  P[0][1]=0;  		P[0][2]=0;                           	   P[0][3]= 0;
	P[1][0]=0;    	  P[1][1]= near/t;  P[1][2]=0;                           	   P[1][3]= 0;
	P[2][0]=0;  	  P[2][1]=0;  		P[2][2]= - (far+near)/(1.0*(far-near));    P[2][3]= -(2.0*far*near)/(far-near);
	P[3][0]=0;  	  P[3][1]=0;  		P[3][2]=-1;                          	   P[3][3]= 0;
}

//===============================Main Function=================================================

int main () {
   init();
   string command;
   double data;
   double t1,t2,t3,t4;
   int numOfTriangles=0;
   // open a file in read mode.
   ifstream infile,infile2,infile3; 
   infile.open("scene.txt"); 

   // open three files in write mode.
   ofstream outfile1,outfile2,outfile3;
   outfile1.open("stage1.txt");
   outfile2.open("stage2.txt");
   outfile3.open("stage3.txt");

   for(int i=0;i<3;i++){	
   		infile >> gLAt.p[i].x >>  gLAt.p[i].y >> gLAt.p[i].z;
   }
   infile >> gLAt.fov >> gLAt.aspectRatio >> gLAt.near >> gLAt.far;




   //=========================Modeling Transformation==================================
   infile >> command ;
   while(command != "end"){
   			if(command ==  "triangle" ){
   				numOfTriangles++;
   				for(int i=0; i<3; i++){
   					for(int j=0; j<3;j++){
   						infile >> triangle[j][i] ;
   					}
					triangle[3][i]=1;
   				}
				handleTriangle();
   				for (int i = 0; i < 3; i++)
			    {
			        for (int j = 0; j < 3; j++)
			        {
		        		outfile1 << fixed << setprecision(7);
		                outfile1 <<  triangle[j][i] << " " ;
			        }
		            outfile1 << endl ;
			    }
			    outfile1 << endl ;
   			}
   				
			if(command ==  "translate" ){
				infile >> t1 >> t2 >> t3;
				handleTranslation(t1,t2,t3);
			}
   				
			if(command ==  "scale" ){
				infile >> t1 >> t2 >> t3;
				handleScaling(t1,t2,t3);
			}
			if(command ==  "rotate" ){
				infile >> t1 >> t2 >> t3 >> t4 ;
				handleRotation(t1,t2,t3,t4);
			}
   				
			if(command ==  "push" ){
				handlePush();
			}
			if(command ==   "pop" ){
				handlePop();
			}

   		infile >> command;
   }




   //==========================View Transformation=====================================
   infile2.open("stage1.txt"); 
   calculateViewTransformationMatrix();
   for (int i = 0; i < numOfTriangles; i++){

   		for (int j = 0; j < 3; j++)
   		{
   			for (int k = 0; k < 3; k++)
   			{
   				infile2 >> triangle[k][j] ;
			}
			triangle[3][j]=1;
   		}

   		for (int j = 0; j < 4; j++)
	    {
	        for (int k = 0; k< 4; k++)
	        {
	            temp[j][k] = 0;
	            for (int l = 0; l < 4; l++)
	            {
	                temp[j][k] += V[j][l] * triangle[l][k];
	            }
	        }
	    }

	    for (int k = 0; k < 4; k++)
	    {
	        for (int j = 0; j < 4; j++)
	        {
	        	if(temp[3][k] == 1 )
	        		triangle[j][k] = temp[j][k];
	        	else //if(temp[3][k] != 0)
	        		triangle[j][k] = temp[j][k]/(temp[3][k]*1.0);
	        }
	    }

	    for (int j = 0; j < 3; j++)
	    {
	        for (int k = 0; k < 3; k++)
	        {	
        		outfile2 << fixed << setprecision(7);
	        	
                outfile2  << triangle[k][j] << " " ;
	        }
            outfile2 << endl ;
	    }
	    outfile2 << endl ;

   }


   //==========================Projection Transformation===============================
   infile3.open("stage2.txt");
   calculateProjectionTransformationMatrix();
   for (int i = 0; i < numOfTriangles; i++){

   		for (int j = 0; j < 3; j++)
   		{
   			for (int k = 0; k < 3; k++)
   			{
   				infile3 >> triangle[k][j] ;
			}
			triangle[3][j]=1;
   		}

   		for (int j = 0; j < 4; j++)
	    {
	        for (int k = 0; k< 4; k++)
	        {
	            temp[j][k] = 0;
	            for (int l = 0; l < 4; l++)
	            {
	                temp[j][k] += P[j][l] * triangle[l][k];
	            }
	        }
	    }

	    for (int k = 0; k < 4; k++)
	    {
	        for (int j = 0; j < 4; j++)
	        {
	        	if(temp[3][k] == 1 )
	        		triangle[j][k] = temp[j][k];
	        	else //if(temp[3][k] != 0)
	        		triangle[j][k] = temp[j][k]/(temp[3][k]*1.0);
	        }
	    }

	    for (int j = 0; j < 3; j++)
	    {
	        for (int k = 0; k < 3; k++)
	        {
        		outfile3 << fixed << setprecision(7);
                outfile3 <<  triangle[k][j] << " " ;
	        }
            outfile3 << endl ;
	    }
	    outfile3 << endl ;

   }

   // close the opened files.
   infile.close();
   infile2.close();
   infile3.close();
   outfile1.close();
   outfile2.close();
   outfile3.close();

   return 0;
}
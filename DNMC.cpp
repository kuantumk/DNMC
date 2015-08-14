// A Monte Carlo simulation code developed 
// to solve the atom diffusion passivation problem. 
// Code developed 2013-2014
// paper published based on this code:
// Applied Physics A
// August 2015, Volume 120, Issue 2, pp 635-639

#include <stdafx.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_N (100000)
#define N_doping (1.0) // doping level, unit: 1E18
#define N_Si ((long)(N_N/100*N_doping))
#define PI (3.141592654)


typedef struct {
	long x, y;
}PosType;


double round(double d) {
    return floor(d + 0.5);
}


double Schrage_random(long &seed) {
	//Schrage random number generator
	long m, a, q, r, z=seed;
	double I=seed, x;
	
	m=2147483647; a=16807; q=127773; r=2836;	//pick a Schrage value
	I=z=a*(z%q)-r*(long)floor(I/q);
	if (I<0) {
		I+=m;
		z+=m;
	}
	x=I/m;
	seed=z;
	return x;
}

double gaussrand() {
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;

	return Z;
}

bool pos_check(long i, long j, PosType *pos, long N) {
	//check mesh position availability: if empty then return true
	int k;
	
	for (k=0; k<N; k++) {
		if (pos[k].x==i && pos[k].y==j) return false;
	}
	return true;
}

	
void assign_pos(PosType *N_pos, PosType *Si_pos, long N_mesh, long &seed) {
	//randomly assign positions for N atoms
	int k, flag;
	long i, j;

	for (k=0; k<N_N; k++) {
		flag=1;
		do {
			i=(long)(Schrage_random(seed)*(double)N_mesh);
			j=(long)(Schrage_random(seed)*(double)N_mesh);
			if (pos_check(i,j,N_pos,k)) {
	  			N_pos[k].x=i;
	  			N_pos[k].y=j;
				flag=0;
			}
		} while (flag==1);
	}

	for (k=0; k<N_Si; k++) {
		flag=1;
		do {
			i=(long)(Schrage_random(seed)*(double)N_mesh);
			j=(long)(Schrage_random(seed)*(double)N_mesh);
			if (pos_check(i,j,N_pos,N_N) && pos_check(i,j,Si_pos,k)) {
				Si_pos[k].x=i;
	  			Si_pos[k].y=j;
				flag=0;
			}
		} while (flag==1);
	}
}


bool neibor_check(long i, long j, PosType *pos, long N, long N_mesh, long &s, long &t) {
	//check if there is an N neibor. Return N atom position (s,t) for N deletion
	int l, r, u, d;
	
	if (i==0) {
		l=N_mesh-1;
		r=i+1;
	}
	else if (i==N-1) {
		l=i-1;
		r=0;
	}
	else {
		l=i-1;
		r=i+1;
	}
	
	if (j==0) {
		u=N_mesh-1;
		d=j+1;
	}
	else if (j==N-1) {
		u=j-1;
		d=0;
	}
	else {
		u=j-1;
		d=j+1;
	}
	
	if (!pos_check(l,u,pos,N)) {
		s=l;
		t=u;
		return true;
	}
	else if (!pos_check(i,u,pos,N)) {
		s=i;
		t=u;
		return true;
	}
	else if (!pos_check(r,u,pos,N)) {
		s=r;
		t=u;
		return true;
	}
	else if (!pos_check(l,j,pos,N)) {
		s=l;
		t=j;
		return true;
	}
	else if (!pos_check(r,j,pos,N)) {
		s=r;
		t=j;
		return true;
	}
	else if (!pos_check(l,d,pos,N)) {
		s=l;
		t=d;
		return true;
	}
	else if (!pos_check(i,d,pos,N)) {
		s=i;
		t=d;
		return true;
	}
	else if (!pos_check(r,d,pos,N)) {
		s=r;
		t=d;
		return true;
	}
	else {
		s=-1;
		t=-1;
		return false;
	}
}


void del_N(long i, long j, PosType *pos, long N) {
	//if N-Si defect complex formed, delete the N position
	int k, l;
	
	for (l=-1, k=0; k<N; k++) {
		if (pos[k].x==i && pos[k].y==j) {
			l=k;
			break;
		}
	}
	if (l!=-1) {
		for (k=l; k<N-1; k++) {
			pos[k].x=pos[k+1].x;
			pos[k].y=pos[k+1].y;
		}
	}
}


double NMobile(double SMDL) {
	long	seed, N_mesh, actual_Si, NDComplex, r, r_cum, x, y, nx, ny, dx, dy, s, t;
	PosType N_pos[N_N], Si_pos[N_Si];
	int		i, j, k, encounter;
	double	u, N_content, ratio, lattice, theta, deg, q_2, q;

	lattice=5.65325;
	N_content=0.02; //N composition=1%
	N_mesh=(long)sqrt(N_N/2/N_content);  //mesh size: N_mesh
	actual_Si=N_Si;
	seed=(long)time(NULL);	// serves as the initial value for random seed
	
	assign_pos(N_pos,Si_pos,N_mesh,seed);

	for (k=0, NDComplex=0; k<N_N; k++) {
		//printf("\nSi number=%i",k);
		x=N_pos[k].x;
		y=N_pos[k].y;
		encounter=0;

		if (neibor_check(x,y,Si_pos,actual_Si,N_mesh,s,t)) {
			//check if Si is already neighbor with N
			NDComplex++;
			del_N(s,t,Si_pos,actual_Si);
			actual_Si--;
			continue;
		}
		
		r_cum=0;
		do {//begin random walk
			u=Schrage_random(seed); //for walking distance
			r=(long)ceil(2*u*SMDL/lattice);
			r_cum+=r;
		
			u=Schrage_random(seed); //for walking angle
			theta=2*PI*u;

			for (j=1; j<=r; j++) {
				dx=(long)round(abs(j*sin(theta)));
				dy=(long)round(abs(j*cos(theta)));

				if (theta<PI/2 || theta>3*PI/2) {
					nx=x+dx;
					if (nx>=N_mesh) nx=nx-N_mesh;
				}
				else {
					nx=x-dx;
					if (nx<0) nx=nx+N_mesh;
				}
				if (theta<PI) {
					ny=y+dy;
					if (ny>=N_mesh) ny=ny-N_mesh;
				}
				else {
					ny=y-dy;
					if (ny<0) ny=ny+N_mesh;
				}
				
				if (neibor_check(nx,ny,Si_pos,actual_Si,N_mesh,s,t)) {
					NDComplex++;
					del_N(s,t,Si_pos,actual_Si);
					actual_Si--;	
					encounter=1;
					break;
				}

				q_2=double((nx-x)^2+(ny-y)^2);
				q=sqrt(q_2);
				if (q >= SMDL) break;
			}

		}while (q < SMDL && !encounter && r_cum<=5000*SMDL);
	}
	ratio=((double)NDComplex)/N_Si;
	return 1-ratio;
}

void main() {
	int i, j, ave_N;
	double SMDL, ratio;
	FILE	*fp;
	
	fp=fopen("Nmobile-2.0.txt", "w+");
	ave_N=10; //average size
	SMDL=5.0;

	do {
		ratio=0;
		for (i=0; i<ave_N; i++) ratio+=NMobile(SMDL);
		ratio=ratio/ave_N;
		fprintf(fp,"%f\t%f\n",SMDL/10,ratio);
		printf("SMDL=%f\n", SMDL);
		SMDL+=5.0;
	} while (SMDL<=500);
	fclose(fp);
}

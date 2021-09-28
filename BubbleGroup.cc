/* 
Program done by Sergio Alonso  
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#define pi 3.141592653589793
using namespace std;

double  Romax =       2.000; 
double  Romin =       0.0;
double  Rspread =     10.000; 
double  Pspread =     0.0; 
double  mult  =       1.0;
int N_Realizations =  14; //  Correponds to the number of classes per school

#define  Nt        20    //   Number of students per class
#define  N_Schools 5131  //   Number of Schools 
#define  IA        250.  //  
#define  t1        100   //  Time to be contagious 
#define  t2        700   //  Time to be less contagious 
#define  t4        1900  // Time to be contagious 
#define  t3        500   //  Time to dtect symptoms -> go home
#define  time_max  50    //   Total time of the simulation in days 
#define  itime_max 5000  //  Total time of the simulation  
#define  tq        1000  //  Time of Quarantine  
#define  NRR       20    //  Number of realizations 
#define Nsymp      500   //  Time to detect because symptoms 
#define Ntrac      1000  //  Time to detect kid due to contact tracing
double dt=1.0* time_max/itime_max;

int Rdistri[Nt];
int Rdistotal=0;

ofstream file_16 ("SCHOOLS.dat");
ofstream file_166 ("ALL_SCHOOLS.dat");
ofstream file_161 ("SCHOOLS_xmgr.dat");
ofstream file_162 ("INFETIONS_SCHOOLS_xmgr.dat");
ofstream file_163 ("INFETIONS_class.dat");
ofstream file_11 ("N_AVEt.dat");
ofstream file_12 ("New_AVE_t.dat");
ofstream file_1 ("N_t.dat");
ofstream file_17 ("N_detection_t.dat");
ofstream file_5 ("His.dat");
ofstream file_gom ("Gom.dat");
ofstream file_numqua ("NumQua_time.dat");
ofstream fileRqu ("Cont_Qua.dat");
ofstream file_61 ("number_infected.dat");
ofstream file_Rdistribucio ("R_distribucio.dat");
ofstream file_RdistribucioFrac ("R_Fraction.dat");


int People[Nt];
int previous[Nt];
int next_inf[Nt];
double RR[Nt];
int Casos[itime_max];
int New_Casos[itime_max];
int Casos_AVE[itime_max];
int New_Casos_AVE[itime_max];
int histo[Nt];
int No=0,N1=0;
int Noo=0;



float f_Gompertz(int time)
{
float Gompertz;
float K=1;  	  
float No=0.001;  
float a=0.005;  
float amuo=log(K/No);
Gompertz= a * K * exp ( -amuo  * exp(-a*time) ) * amuo * exp(-a*time) ;
return Gompertz;
}




int main()
{
  double t=0;
  int it=0,jt=0,kt=0,kkt=0;
  double random=0.0;
  double random2=0.0; 
  int k=0;
  int quarantine=0;
  int timeini;
  int N2=0,N=0;
  int N2tau=0,Ntau=0;
  int ischool=0;
  double p1=0.0;
  srand (time(NULL));
  int irr=0;
  int ohhh=0;
  int classfull=0;

  int num_Qua[itime_max+tq];
   for (k=0;k<itime_max+tq;k++) {num_Qua[k]=0;}
   for (k=0;k<Nt;k++) {Rdistri[k]=0;}


  for (irr=0;irr<NRR;irr++) {

// For all school
  int Contagious_Initial_tot=0;
  int Contagious_outside_tot=0;
  int Contagious_inside_tot=0;
  int Contagious_family_tot=0;
  int Contagious_assimp_tot=0;
  int Contagious_simp_tot=0;
  int Contagious_tracing_tot=0;
  int nq=0;


  for (ischool=0;ischool<N_Schools;ischool++) {
  for (k=0;k<Nt;k++) {histo[k]=0;}

// For each school
  int Contagious_Initial=0;
  int Contagious_outside=0;
  int Contagious_inside=0;
  int Contagious_family=0;
  int Contagious_assimp=0;
  int Contagious_simp=0;
  int Contagious_tracing=0;
  int Contagious_nodetect=0;


 for (it=0;it<itime_max;it++) {Casos_AVE[it]=0;New_Casos_AVE[it]=0;}

  for (jt=0;jt<N_Realizations;jt++) {
  quarantine=0;
  

//Inititial Condition
  Noo=0;

  random = rand()/(RAND_MAX+1.0);
  if (random< 0.3* IA*Nt/100000) { Noo++;
  			  random = rand()/(RAND_MAX+1.0);
  			  if (random< 0.3 *IA*Nt/100000) { Noo++;
				random = rand()/(RAND_MAX+1.0);
  			  	if (random< IA*Nt/100000) { Noo++;}
				}
			}

  Contagious_Initial+=Noo;

//generation of the initial contagius children
  for (k=0;k<Noo;k++) {random = rand()/(RAND_MAX+1.0); timeini=(random*t2); People[k]=timeini; previous[k]=-1; next_inf[k]=0; random = rand()/(RAND_MAX+1.0); RR[k]=(Romin+random*(Romax-Romin));
  random = rand()/(RAND_MAX+1.0); if (random<Pspread) {RR[k]=Rspread;}    }   

  N=Noo; Ntau=0;
  t=0;


//time evolution
  for (it=0;it<itime_max;it=it+1) {
  t=t+dt;
  N2=N;
  N2tau=Ntau;

  if (quarantine==tq) {quarantine=0;}
  
  if  ( quarantine>=1 && quarantine<tq) {quarantine+=1; for (k=0;k<N;k++) {People[k]=People[k]+1;} num_Qua[it]+=1;}
  else {
  for (k=0;k<N;k++) {
  People[k]=People[k]+1;
   
  			
  random = rand()/(RAND_MAX+1.0);
  random2 = rand()/(RAND_MAX+1.0);
  p1 = 1.0 * f_Gompertz (People[k]-t1);


	if (People[k]==Nsymp ) { 	  
		if (random < 0.3) { ohhh=0; for (kkt=0;kkt<N;kkt++) {if (People[kkt]<t4) {ohhh++; 
		}} fileRqu<<nq<<" "<<ohhh<< " " << N << " " << jt << " " << ischool<<"\n";	Rdistri[ohhh]++;
				    N2tau++;    
				quarantine=1; nq++; 
				Contagious_simp+=1; 
				} //kid is detected after 5 days because of simptoms! or because contact tracing outside school !
	else {Contagious_assimp+=1; } //kid is not detected after 5 days because is assimptomatic!
	}
	
	if (previous[k]<=-1) { 
	if (People[k]==Nsymp ) {	  
		if (random < 0.7 && random > 0.3) { ohhh=0; for (kkt=0;kkt<N;kkt++) {if (People[kkt]<t4) {ohhh++;
		}} fileRqu<<nq<<" "<<ohhh<< " " << N << " " << jt << " " << ischool<<"\n"; Rdistri[ohhh]++;	
				    N2tau++;    
				quarantine=1; nq++; 
				Contagious_tracing+=1; } //kid is detected after 5 days because of simptoms! or because contact tracing outside school !
	else { Contagious_nodetect+=1; } //kid is not detected after 5 days because is assimptomatic!
	}
	}	
	else {cout << previous[k] << "\n";
	if (People[k]==Ntrac ) { 	  
		if (random < 0.57) { ohhh=0; for (kkt=0;kkt<N;kkt++) {if (People[kkt]<t4) {ohhh++;
		}} fileRqu<<nq<<" "<<ohhh<< " " << N << " " << jt << " " << ischool<<"\n";	Rdistri[ohhh]++;	
				    N2tau++;    
				quarantine=1; nq++; 
				Contagious_tracing+=1; } //kid is detected after 5 days because of simptoms! or because contact tracing outside school !
	else { Contagious_nodetect+=1; } //kid is not detected after 5 days because is assimptomatic!
	}
	}
	
	
  if (People[k]<t4) { 


	  if (People[k]>t1) {  
		random = rand()/(RAND_MAX+1.0);          //Kid infects some familiar outside school
		if (random < (0.70  )) {
			p1 = 1.0 * f_Gompertz (People[k]-t1);   //R outside = 1
		  	random = rand()/(RAND_MAX+1.0);
		  		if (random < p1) {
				Contagious_family+=1;  }    
			}
		else{  
			p1 = RR[k] * f_Gompertz (People[k]-t1);
		  	random = rand()/(RAND_MAX+1.0);
		  	random2 = rand()/(RAND_MAX+1.0);
		  	if (random < p1  && random2 < 1.0*(Nt-N)/Nt ) {
				People[N2]=1; 
				previous[N2]=k; 
				next_inf[k]+=1; 
				random = rand()/(RAND_MAX+1.0); 
				RR[N2]=(Romin+random*(Romax-Romin));   random = rand()/(RAND_MAX+1.0); if (random<Pspread) {RR[k]=Rspread;} 
				histo[k]++;    
				N2=N2+1; 		
				Contagious_inside+=1;}    //kid infects other kid of the school 34 %
		}

}
	  }
   }
   }
 if  ( quarantine==0) {
//kid can be infected outside the school
random = rand()/(RAND_MAX+1.0);
random2 = rand()/(RAND_MAX+1.0);		  	
if (random < 0.7 * IA*Nt/(100000*1400)  && random2 < 1.0*(Nt-N)/Nt ) { People[N2]=1; previous[N2]=-2; next_inf[k]=0; random = rand()/(RAND_MAX+1.0); RR[N2]=(Romin+random*(Romax-Romin));   
					random = rand()/(RAND_MAX+1.0); if (random<Pspread) {RR[k]=Rspread;} 
				        N2=N2+1;Contagious_outside+=1;}  

}
		

  Casos[it]=N2tau-Ntau;
  Casos_AVE[it]+=Ntau;
  New_Casos_AVE[it]+=N2tau-Ntau;
  N=N2; Ntau=N2tau;
  if (N>=Nt) {classfull++;}

}


for (k=0;k<N;k++) {file_61 << (ischool+1)*10000+ jt*1000 +  k << " " <<  next_inf[k]   << "\n";}
for (k=0;k<N;k++) {if (histo[k]>=0)  file_5 <<  histo[k] << "\n";}



file_163  << N << "\n";

}
file_16 << ischool  << ", " <<  Contagious_Initial << ", " <<  Contagious_outside   << ", " <<  Contagious_inside  << ", " <<  Contagious_family << ", " << Contagious_simp << ", " <<  Contagious_assimp << ", " <<  Contagious_tracing <<      "\n";

file_161 << ischool  << " " <<  Contagious_Initial << " " <<  Contagious_outside   << " " <<  Contagious_inside  << " " <<  Contagious_family << " " << Contagious_simp << " " <<  Contagious_assimp  << " " <<  Contagious_tracing  << "\n";

file_162 << ischool  << " " <<  Contagious_Initial + Contagious_outside   +  Contagious_inside << " " << Contagious_simp + Contagious_assimp << "\n";


Contagious_Initial_tot += Contagious_Initial;
Contagious_outside_tot += Contagious_outside;
Contagious_inside_tot  += Contagious_inside;
Contagious_family_tot  += Contagious_family;
Contagious_assimp_tot  += Contagious_assimp;
Contagious_simp_tot    += Contagious_simp;
Contagious_tracing_tot += Contagious_tracing;


}

file_166 << ischool  << " " <<  Contagious_Initial_tot << " " <<  Contagious_outside_tot   << " " <<  Contagious_inside_tot  << " " <<  Contagious_family_tot << " " << Contagious_simp_tot << " " <<  Contagious_assimp_tot << " " <<  Contagious_tracing_tot << " " << nq <<     "\n";


}

for (kkt=0;kkt<itime_max;kkt++) {file_numqua  << kkt << " " << num_Qua[kkt] /(1.0 * NRR)<< "\n";}
for (it=0;it<1900;it+=10) { file_gom << it+t1  <<  " " <<  f_Gompertz (it) << "\n";  }
for (it=1;it<Nt;it++) { Rdistotal+= Rdistri[it]; }
for (it=1;it<Nt;it++) { file_Rdistribucio << it  <<  " " <<  Rdistri[it]  << "\n"; file_RdistribucioFrac << it  <<  " " <<  1.0*Rdistri[it] / Rdistotal << "\n"; }

file_1.close(); 
file_5.close(); 

return 0;  
} 






















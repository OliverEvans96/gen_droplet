

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TMath.h"

using namespace std;

int verbose = 0;
int fWrite = 1;
int fShift = 0;
int fShiftCentTozero = 1;
int fFullSphere = 0;
int fFullCylinder = 0;
int fSlab = 0;
int fCombine = 0;
double Rmax = 50;


//For PMMA (Mesfin) 
// double xSlabLo = 2.0;
// double xSlabHi = 73.0;
// double ySlabLo = 0.0;
// double ySlabHi = 75.0;
// double zSlabLo = 0.0;
// double zSlabHi = 25.0;
 
// double xSlabLo = 2.5;
// double xSlabHi = 73.0;
// double ySlabLo = -0.5;
// double ySlabHi = 75.5;
// double zSlabLo = 0.0;
// double zSlabHi = 25.0;

// double xSlabLo = 0.0;
// double xSlabHi = 75.0;
// double ySlabLo = 0.0;
// double ySlabHi = 75.0;
// double zSlabLo = 0.0;
// double zSlabHi = 25.0;

double xSlabLo = 0.0;
double xSlabHi = 36.36;
double ySlabLo = 0.0;
double ySlabHi = 36.36;
double zSlabLo = 0.0;
double zSlabHi = 27.52;

// //For Sapphire
// double xSlabLo = -0.801891;	
// double xSlabHi = 53.1289;
// double ySlabLo = -0.75391;	
// double ySlabHi = 50.1901;
// double zSlabLo = -50.0;
// double zSlabHi = 50.0;

//For excess chem pot
// double xSlabLo = -15.0;//-1.0099380075100193e+01;//-15;
// double xSlabHi = 15.0;//1.0099380075100193e+01;//15;
// double ySlabLo = -15.0;//-1.2948472242728267e+01;//-15;
// double ySlabHi = 15.0;//1.2948472242728267e+01;//15;
// double zSlabLo = -15.0;//-28.5;//-15;
// double zSlabHi = 15.0;//28.5;//15;
// double yCylinderLo = 0.9;
// double yCylinderHi = 49.1;
//double yCylinderLo = -0.2;//-0.75398; 
//double yCylinderHi = 50.0;
double yCylinderLo = -1.6267196;      
double yCylinderHi = 108.7384636;

int run1 = 20000,run2 = 1000; 
int dump_steps1 = 1000, dump_steps2 = 100;
int n_time_steps1, n_time_steps2;

int time_step;

double M;
int N_mol1,N_mol2,N_mol_sp;  
int N_atoms_mol1 = 3, N_atoms_mol2 = 3;
int N_atoms1,N_atoms2;
int Nbonds_mol1 = 2,Nbonds_mol2 = 2;
int Nangles_mol1 = 1,Nangles_mol2 = 1;

int N_bonds1,N_bonds2;
int N_angles1,N_angles2;
int N_dihedrals1, N_dihedrals2;
int N_impropers1, N_impropers2;

int N_atom_types1 = 2, N_atom_types2 = 2;
int N_bond_types1 = 1, N_bond_types2 = 1;
int N_angle_types1 = 1, N_angle_types2 = 1;

int N_atoms_sp;
int N_bonds_sp;
int N_angles_sp;
int N_dihedrals_sp = 0;
int N_impropers_sp = 0;

int N_atoms;
int id, type, mol;
double xs, ys, zs; //scaled position with in box
double ix, iy, iz; //number of times boundary is crossed
double xu, yu, zu; //actual position
double x, y, z;    //position wrapped into box


int Check = 0;
int CheckDroplet = 0;
int CheckSlab = 0;
int UseActualPos = 0;

//box sides

double xlo,xhi,Lx;
double ylo,yhi,Ly;  
double zlo,zhi,Lz;
double volume;

//1st box sides
double xb1lo,xb1hi;
double yb1lo,yb1hi;  
double zb1lo,zb1hi;
//center of 1st box
double xb1c, yb1c, zb1c;

//2nd box sides
double xb2lo,xb2hi;
double yb2lo,yb2hi;  
double zb2lo,zb2hi;
//center of 2nd box
double xb2c, yb2c, zb2c;

double dtz = 0.0;  //dtz (in Angstrom) allows for z separation between upper and lower box


//  Molecule Water;

int* Time_step1;
int ** Atom_id1,** Atom_type1,** Atom_mol1;
double** Atom_x1,** Atom_y1,** Atom_z1;
double** Mol_xcm1,** Mol_ycm1,** Mol_zcm1;
double* Vol1,*rho1, *N_den1;

int* Time_step2;
int ** Atom_id2,** Atom_type2,** Atom_mol2;
double** Atom_x2,** Atom_y2,** Atom_z2;
double** Mol_xcm2,** Mol_ycm2,** Mol_zcm2;
double* Vol2,*rho2, *N_den2;


int* Atom_spid;
int* Atom_spmol;
int* Atom_sptype;
double* Atom_spx;   
double* Atom_spy;
double* Atom_spz;

int* d_id;
int* d_type;
int* d_mol;
   
double* d_xs;
double* d_ys;
double* d_zs;

//double charge[2] = {-0.798, 0.399}; // PMMA for Mesfin
double charge[2] = {-0.8476,0.4238};  //SPCE charge[0]=Oxygen, charge[1]=Hydrogen
double Mass[2] = {15.9994,1.00794};   //Mass[0]=Oxygen, Mass[1]=Hydrogen
double rho_exact = 1.0e-24; //density of water (in gm/A^3)

double NA = 6.023e+23; //Avogadro's number

// struct Molecule
// {
//   int N_atoms;
//   int atom_id[3];
//   int atom_type[3];
//   int mol_id;
//   double x[3], y[3], z[3];
//   double vx[3],vy[3],vz[3];

// };

void ShiftIntoBox(int t_step){

  //Symmetry is broken in x, y, z, i.e, no folding back of atoms
  // that would have been out side of simulation box

  double Lx = (xhi-xlo);
  double Ly = (yhi-ylo);
  double Lz = (zhi-zlo);
  int id  = 0;

  for(int ml = 0; ml < N_mol1; ml++){
    Atom_x1[t_step][id] = Atom_x1[t_step][id] - Lx*((int)(Mol_xcm1[t_step][ml]/Lx));
    Atom_y1[t_step][id] = Atom_y1[t_step][id] - Ly*((int)(Mol_ycm1[t_step][ml]/Ly));
    Atom_z1[t_step][id] = Atom_z1[t_step][id] - Lz*((int)(Mol_zcm1[t_step][ml]/Lz));
    id++;
  }
}



void CalcBoxCenter(){

  cout<<"CalcBoxCenter()"<<endl;
  //center of 1st box

  xb1c = 0.5*(xb1lo + xb1hi);
  yb1c = 0.5*(yb1lo + yb1hi);
  zb1c = 0.5*(zb1lo + zb1hi);
  cout<<"xb1c = "<<xb1c<<" , yb1c = "<<yb1c<<" , zb1c = "<<zb1c<<endl;

  //center of 2nd box
  xb2c = 0.5*(xb2lo + xb2hi);
  yb2c = 0.5*(yb2lo + yb2hi);
  zb2c = 0.5*(zb2lo + zb2hi);
  cout<<"xb2c = "<<xb2c<<" , yb2c = "<<yb2c<<" , zb2c = "<<zb2c<<endl;

  cout<<"Done in CalcBoxCenter()"<<endl;
}


void ShiftCenterToZero(){

 //geometric box center
 double xbc = 0.5*(xb1lo + xb1hi);
 double ybc = 0.5*(yb1lo + yb1hi);
 double zbc = 0.5*(zb1lo + zb1hi);
 
  N_atoms_sp = N_mol_sp*N_atoms_mol1;

  for(int m = 0; m < N_atoms_sp; m++){
    Atom_spx[m] = Atom_spx[m]-xbc;
    Atom_spy[m] = Atom_spy[m]-ybc;
    if(fSlab)
      Atom_spz[m] = Atom_spz[m]-zbc;
  }

 xb1lo = xb1lo - xbc;
 yb1lo = yb1lo - ybc;
 if(fSlab)
   zb1lo = zb1lo - zbc;
 xb1hi = xb1hi - xbc;
 yb1hi = yb1hi - ybc;
 if(fSlab)
   zb1hi = zb1hi - zbc;
}

void DuplicateData(){

  cout<<"********************"<<endl;
  cout<<"Duplicating Data()"<<endl;
  cout<<"********************"<<endl;

 for(int i = 0; i < 1; i++){
   //Fill original box
   for(int m = 0; m < N_atoms1; m++){

     double x = Atom_x1[i][m];
     double y = Atom_y1[i][m];
     double z = Atom_z1[i][m];
     
     if(m%3 == 0){
       double xH1 = Atom_x1[i][m+1];
       double xH2 = Atom_x1[i][m+2];
       if( fabs(x - xH1) > 10. && x - xH1 > 0)
	 Atom_x1[i][m+1] = Atom_x1[i][m+1] + Lx;
       if( fabs(x - xH1) > 10. && x - xH1 < 0)
	 Atom_x1[i][m+1] = Atom_x1[i][m+1] - Lx;
       
       if( fabs(x - xH2) > 10. && x - xH2 > 0)
	 Atom_x1[i][m+2] = Atom_x1[i][m+2] + Lx;
       if( fabs(x - xH2) > 10. && x - xH2 < 0)
	 Atom_x1[i][m+2] = Atom_x1[i][m+2] - Lx;
       
       double yH1 = Atom_y1[i][m+1];
       double yH2 = Atom_y1[i][m+2];
       if( fabs(y - yH1) > 10. && y - yH1 > 0)
	 Atom_y1[i][m+1] = Atom_y1[i][m+1] + Ly;
       if( fabs(y - yH1) > 10. && y - yH1 < 0)
	 Atom_y1[i][m+1] = Atom_y1[i][m+1] - Ly;
       
       if( fabs(y - yH2) > 10. && y - yH2 > 0)
	 Atom_y1[i][m+2] = Atom_y1[i][m+2] + Ly;
       if( fabs(y - yH2) > 10. && y - yH2 < 0)
	 Atom_y1[i][m+2] = Atom_y1[i][m+2] - Ly;
       
       double zH1 = Atom_z1[i][m+1];
       double zH2 = Atom_z1[i][m+2];
       if( fabs(z - zH1) > 10. && z - zH1 > 0)
	 Atom_z1[i][m+1] = Atom_z1[i][m+1] + Lz;
       if( fabs(z - zH1) > 10. && z - zH1 < 0)
	 Atom_z1[i][m+1] = Atom_z1[i][m+1] - Lz;
       
       if( fabs(z - zH2) > 10. && z - zH2 > 0)
	 Atom_z1[i][m+2] = Atom_z1[i][m+2] + Lz;
       if( fabs(z - zH2) > 10. && z - zH2 < 0)
	 Atom_z1[i][m+2] = Atom_z1[i][m+2] - Lz;
     }
     
   }
 }

 //geometric box center
 double xbc = 0.5*(xb1lo + xb1hi);
 double ybc = 0.5*(yb1lo + yb1hi);
 double zbc = 0.5*(zb1lo + zb1hi);

 int id = 0;
 
 int N = 8*N_atoms1;

 int* d_id = new int[N];
 int* d_type = new int[N];
 int* d_mol = new int[N];
   
 double* d_xs = new double[N];
 double* d_ys = new double[N];
 double* d_zs = new double[N];
 
 for(int i = 0; i < 1; i++){
   //Fill original box
   for(int m = 0; m < N_atoms1; m++){
     
     d_id[id] = Atom_id1[i][m];
     d_type[id] = Atom_type1[i][m];
     d_mol[id] = Atom_mol1[i][m];
     
     d_xs[id] = Atom_x1[i][m];
     d_ys[id] = Atom_y1[i][m];
     d_zs[id] = Atom_z1[i][m];
     
     id++;
     
   }
   //Duplicate original box in the x-direction
   for(int m = 0; m < N_atoms1; m++){
     
     d_id[id] = Atom_id1[i][m];
     d_type[id] = Atom_type1[i][m];
     d_mol[id] = Atom_mol1[i][m];
     
     d_xs[id] = Atom_x1[i][m] + (xb1hi - xb1lo);
     d_ys[id] = Atom_y1[i][m];
     d_zs[id] = Atom_z1[i][m];
      
      id++;
      
    }
   //duplicate original box in the y direction
    for(int m = 0; m < N_atoms1; m++){
   
      d_id[id] = Atom_id1[i][m];
      d_type[id] = Atom_type1[i][m];
      d_mol[id] = Atom_mol1[i][m];
      
      d_xs[id] = Atom_x1[i][m];
      d_ys[id] = Atom_y1[i][m] + (yb1hi - yb1lo);
      d_zs[id] = Atom_z1[i][m];
      
      id++;
      
    } 
    //duplicate the original box in the x and y direction
    for(int m = 0; m < N_atoms1; m++){
   
      d_id[id] = Atom_id1[i][m];
      d_type[id] = Atom_type1[i][m];
      d_mol[id] = Atom_mol1[i][m];
      
      d_xs[id] = Atom_x1[i][m] + (xb1hi - xb1lo);
      d_ys[id] = Atom_y1[i][m] + (yb1hi - yb1lo);
      d_zs[id] = Atom_z1[i][m];
      
      id++;
      
    }
    //duplicate original box in z-direction
    for(int m = 0; m < N_atoms1; m++){
   
      d_id[id] = Atom_id1[i][m];
      d_type[id] = Atom_type1[i][m];
      d_mol[id] = Atom_mol1[i][m];
      
      d_xs[id] = Atom_x1[i][m];
      d_ys[id] = Atom_y1[i][m];
      d_zs[id] = Atom_z1[i][m] + (zb1hi - zb1lo);
      
      id++;
      
    } 
     //duplicate original box in x and z-directions
    for(int m = 0; m < N_atoms1; m++){
   
      d_id[id] = Atom_id1[i][m];
      d_type[id] = Atom_type1[i][m];
      d_mol[id] = Atom_mol1[i][m];
      
      d_xs[id] = Atom_x1[i][m] + (xb1hi - xb1lo);
      d_ys[id] = Atom_y1[i][m];
      d_zs[id] = Atom_z1[i][m] + (zb1hi - zb1lo);
      
      id++;
      
    }
    //duplicate original box in y and z-direction
    for(int m = 0; m < N_atoms1; m++){
   
      d_id[id] = Atom_id1[i][m];
      d_type[id] = Atom_type1[i][m];
      d_mol[id] = Atom_mol1[i][m];
      
      d_xs[id] = Atom_x1[i][m];
      d_ys[id] = Atom_y1[i][m] + (yb1hi - yb1lo);
      d_zs[id] = Atom_z1[i][m] + (zb1hi - zb1lo);
      
      id++;
      
    }
     //duplicate original box in x,y and z-direction
    for(int m = 0; m < N_atoms1; m++){
   
      d_id[id] = Atom_id1[i][m];
      d_type[id] = Atom_type1[i][m];
      d_mol[id] = Atom_mol1[i][m];
      
      d_xs[id] = Atom_x1[i][m] + (xb1hi - xb1lo);
      d_ys[id] = Atom_y1[i][m] + (yb1hi - yb1lo);
      d_zs[id] = Atom_z1[i][m] + (zb1hi - zb1lo);
      
      id++;
      
    }
    
  }


 cout<<"center of original box"<<endl;
 cout<<"xb1c = "<<xb1c<<" , yb1c = "<<yb1c<<" , zb1c = "<<zb1c<<endl;
 xb1lo = xb1lo;
 xb1hi = xb1hi + (xb1hi - xb1lo);
 yb1lo = yb1lo;
 yb1hi = yb1hi + (yb1hi - yb1lo);
 zb1lo = zb1lo;
 zb1hi = zb1hi + (zb1hi - zb1lo);


 Lx = xb1hi-xb1lo;
 Ly = yb1hi-yb1lo;
 Lz = zb1hi-zb1lo;	

 //center of original box

 xb1c = 0.5*(xb1lo + xb1hi);
 yb1c = 0.5*(yb1lo + yb1hi);
 zb1c = 0.5*(zb1lo + zb1hi);
 cout<<"center of duplicated box"<<endl;
 cout<<"xb1c = "<<xb1c<<" , yb1c = "<<yb1c<<" , zb1c = "<<zb1c<<endl;
 
 //shift box so that center is at (0,0,0)
 xb2lo = xb1lo - xb1c;
 xb2hi = xb1hi - xb1c;
 yb2lo = yb1lo - yb1c;
 yb2hi = yb1hi - yb1c;
 zb2lo = zb1lo - zb1c;
 zb2hi = zb1hi - zb1c;
 
 //center of shifted box
 xb2c = 0.5*(xb2lo + xb2hi);
 yb2c = 0.5*(yb2lo + yb2hi);
 zb2c = 0.5*(zb2lo + zb2hi);
 cout<<"center of shifted box"<<endl;
 cout<<"xb2c = "<<xb2c<<" , yb2c = "<<yb2c<<" , zb2c = "<<zb2c<<endl;

 delete [] Atom_id1[0];
 delete [] Atom_type1[0];
 delete [] Atom_mol1[0];
 delete [] Atom_x1[0];
 delete [] Atom_y1[0];
 delete [] Atom_z1[0];

 N_atoms1 = N;


 Atom_id1[0] = new int[N_atoms1];
 Atom_type1[0] = new int[N_atoms1];
 Atom_mol1[0] = new int[N_atoms1];
 Atom_x1[0] = new double[N_atoms1];
 Atom_y1[0] = new double[N_atoms1];
 Atom_z1[0] = new double[N_atoms1];
 //shift atom positions to within shifted box
 for(int m = 0; m < N_atoms1; m++){
   
   Atom_id1[0][m] = d_id[m];
   Atom_type1[0][m] = d_type[m];
   Atom_mol1[0][m] = d_mol[m];
   
   Atom_x1[0][m] = d_xs[m];
   Atom_y1[0][m] = d_ys[m];
   Atom_z1[0][m] = d_zs[m];
 }

 
	
 //  //shift atom positions to within shifted box
//  for(int m = 0; m < N; m++){
//    d_xs[m] = d_xs[m]-xb1c;
//    d_ys[m] = d_ys[m]-yb1c;
//    d_zs[m] = d_zs[m]-zb1c;
//  }

 cout<<"********************"<<endl;
 cout<<"Done Duplicating Data()"<<endl;
 cout<<"********************"<<endl;

}

void WriteDuplicateData(string name, int t_step){


  fstream outfile;

  outfile.open(name.c_str(),ios::out);
 
 int N_atoms_tot = N_atoms1;

 N_mol1 = N_atoms_tot/N_atoms_mol1;

 cout<<"N_mol1 = "<<N_mol1<<endl;

 int N_bonds_tot = 2*N_mol1;
 int N_angles_tot = N_mol1;
 int N_dihedrals_tot = 0;
 int N_impropers_tot = 0;
 
 

 double dx,dy,dz,dr;
 
 int i = t_step, atid;
 int ix = 0,iy = 0,iz = 0;
 
 N_mol_sp = 0;
  
 for(int j = 0; j < N_mol1; j++){
    
   //use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = d_xs[atid-1] - xb2c;
   dy = d_ys[atid-1] - yb2c;
   dz = d_zs[atid-1] - zb2c;
   
   if(fabs(dx) < Rmax && fabs(dy) < Rmax && fabs(dz) < Rmax){
     dr = sqrt(dx*dx + dy*dy + dz*dz);
     
     if(dr < Rmax) N_mol_sp++;
   }
 }
 
 cout<<"N_droplet_mol = "<<N_mol_sp<<endl;
 cout<<"N_droplet_atom = "<<N_mol_sp*N_atoms_mol1<<endl;

 Atom_spid = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_spmol = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_sptype = new int[N_mol_sp*N_atoms_mol1]; 
 
 Atom_spx = new double[N_mol_sp*N_atoms_mol1];   
 Atom_spy = new double[N_mol_sp*N_atoms_mol1];  
 Atom_spz = new double[N_mol_sp*N_atoms_mol1];  
 
 
 int nid = 0,mid = 0,type;
 double ch,x,y,z;
 
 for(int j = 0; j < N_mol1; j++){
   
   //use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   
   dx = d_xs[atid-1] - xb2c;
   dy = d_ys[atid-1] - yb2c;
   dz = d_zs[atid-1] - zb2c;
   
    if(fabs(dx) < Rmax && fabs(dy) < Rmax && fabs(dz) < Rmax){
      dr = sqrt(dx*dx + dy*dy + dz*dz);
      
      if(dr < Rmax){
	
	if(CheckDroplet)
	  cout<<"Mol: "<<j+1<<" , dx = "<<dx<<", dy = "<<dy<<" , dz = "<<dz<<", dr = "<<dr<<" , Rmax = "<<Rmax<<endl;
	
	for(int k = 0; k < 3; k++){
	  
	  id  = N_atoms_mol1*j + k + 1;
	  type = d_type[id-1];
	  mol = d_mol[id-1]; 
	  
	  x = d_xs[id-1];   
	  y = d_ys[id-1];
	  z = d_zs[id-1]; 

	  if(k == 0){
	    double xH1 = d_xs[id];
	    double xH2 = d_xs[id+1];
	    if( fabs(x - xH1) > 10. && x - xH1 > 0)
	      d_xs[id] = d_xs[id] + 0.5*(xb2hi - xb2lo);
	    if( fabs(x - xH1) > 10. && x - xH1 < 0)
	      d_xs[id] = d_xs[id] - 0.5*(xb2hi - xb2lo);
	    
	    if( fabs(x - xH2) > 10. && x - xH2 > 0)
	      d_xs[id+1] = d_xs[id+1] + 0.5*(xb2hi - xb2lo);
	    if( fabs(x - xH2) > 10. && x - xH2 < 0)
	      d_xs[id+1] = d_xs[id+1] - 0.5*(xb2hi - xb2lo);

	    double yH1 = d_ys[id];
	    double yH2 = d_ys[id+1];
	    if( fabs(y - yH1) > 10. && y - yH1 > 0)
	      d_ys[id] = d_ys[id] + 0.5*(yb2hi - yb2lo);
	    if( fabs(y - yH1) > 10. && y - yH1 < 0)
	      d_ys[id] = d_ys[id] - 0.5*(yb2hi - yb2lo);
	    
	    if( fabs(y - yH2) > 10. && y - yH2 > 0)
	      d_ys[id+1] = d_ys[id+1] + 0.5*(yb2hi - yb2lo);
	    if( fabs(y - yH2) > 10. && y - yH2 < 0)
	      d_ys[id+1] = d_ys[id+1] - 0.5*(yb2hi - yb2lo);

	    double zH1 = d_zs[id];
	    double zH2 = d_zs[id+1];
	    if( fabs(z - zH1) > 10. && z - zH1 > 0)
	      d_zs[id] = d_zs[id] + 0.5*(zb2hi - zb2lo);
	    if( fabs(z - zH1) > 10. && z - zH1 < 0)
	      d_zs[id] = d_zs[id] - 0.5*(zb2hi - zb2lo);
	    
	    if( fabs(z - zH2) > 10. && z - zH2 > 0)
	      d_zs[id+1] = d_zs[id+1] + 0.5*(zb2hi - zb2lo);
	    if( fabs(z - zH2) > 10. && z - zH2 < 0)
	      d_zs[id+1] = d_zs[id+1] - 0.5*(zb2hi - zb2lo);
	  }

	  Atom_spid[nid] = nid;
	  Atom_spmol[nid] = j+1;
	  Atom_sptype[nid] = type;
	  Atom_spx[nid] = x;   
	  Atom_spy[nid] = y;
	  Atom_spz[nid] = z;
	  nid++; 
	  if(CheckDroplet)
	    cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	}
      }
    }
    
 }
  
  if(verbose){
    cout<<"LAMMPS Atom File"<<endl;
    cout<<endl;
    cout<<"          "<<N_atoms_tot<<" atoms"<<endl;
    cout<<"          "<<N_bonds_tot<<" bonds"<<endl;
    cout<<"          "<<N_angles_tot<<" angles"<<endl;
    cout<<"          "<<N_dihedrals_tot<<" dihedrals"<<endl;
    cout<<"          "<<N_impropers_tot<<" impropers"<<endl;
    cout<<endl;
    cout<<"          "<<N_atom_types1<<" atom types"<<endl;
    cout<<"          "<<N_bond_types1<<" bond types"<<endl;
    cout<<"          "<<N_angle_types1<<" angle types"<<endl;
    cout<<endl;
    
    cout<<"     "<<xb1lo<<"    "<<xb1hi + (xb1hi - xb1lo)<<"  xlo xhi"<<endl;
    cout<<"     "<<yb1lo<<"    "<<yb1hi + (yb1hi - yb1lo)<<"  ylo yhi"<<endl;
    cout<<"     "<<zb1lo<<"    "<<zb1hi + (zb1hi - zb1lo) <<"  zlo zhi"<<endl;
    
    cout<<endl;  
    cout<<"Masses"<<endl;
    cout<<endl;
    cout<<"\t"<<1<<"\t"<<Mass[0]<<endl;
    cout<<"\t"<<2<<"\t"<<Mass[1]<<endl;
  }

  if(fWrite){
   
    outfile<<"LAMMPS Atom File"<<endl;
    outfile<<endl;
    outfile<<"          "<<N_atoms_tot<<" atoms"<<endl;
    outfile<<"          "<<N_bonds_tot<<" bonds"<<endl;
    outfile<<"          "<<N_angles_tot<<" angles"<<endl;
    outfile<<"          "<<N_dihedrals_tot<<" dihedrals"<<endl;
    outfile<<"          "<<N_impropers_tot<<" impropers"<<endl;
    outfile<<endl;
    outfile<<"          "<<N_atom_types1<<" atom types"<<endl;
    outfile<<"          "<<N_bond_types1<<" bond types"<<endl;
    outfile<<"          "<<N_angle_types1<<" angle types"<<endl;
    outfile<<endl;
    
    outfile<<"     "<<xb2lo<<"    "<<xb2hi<<"  xlo xhi"<<endl;
    outfile<<"     "<<yb2lo<<"    "<<yb2hi<<"  ylo yhi"<<endl;
    outfile<<"     "<<zb2lo<<"    "<<zb2hi + (zb1hi - zb1lo) <<"  zlo zhi"<<endl;
  
    outfile<<endl;  
    outfile<<"Masses"<<endl;
    outfile<<endl;
    outfile<<"\t"<<1<<"\t"<<Mass[0]<<endl;
    outfile<<"\t"<<2<<"\t"<<Mass[1]<<endl;
    
 
  }
  
  //##############################
  //## generate atom list
  //##############################
	  
 
  cout<<endl;
  cout<<"Atoms"<<endl; 
  cout<<endl;
  
  if(fWrite){
    outfile<<endl;
    outfile<<"Atoms"<<endl;
    outfile<<endl;
  }


  //write out droplet first
  for(int m = 0; m < N_atoms_tot; m++){
   
    nid = m+1;
    if(m%N_atoms_mol1 == 0) mid++;
    type = d_type[m];
    x  = d_xs[m];   
    y  = d_ys[m];
    z  = d_zs[m];

    if(m%3 == 0)
      ch = charge[0];
    else
      ch = charge[1];

    if(verbose){
      cout<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	  <<x<<"  "<<y<<"  "<<z<<"  0  0  0"<<endl;
    }
    if(fWrite){
      outfile<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	     <<x<<"  "<<y<<"  "<<z<<"  0  0  0"<<endl;
    }	
  
//       1    1  1 -0.8472   12.12456   28.09298   22.27452  0  1  0
//       2    1  2  0.4236   12.53683   28.75606   22.89928  0  1  0
//       3    1  2  0.4236   11.49482   28.56390   21.65678  0  1  0

  }


// //########################
// //## generate bond list
// //########################
  

  cout<<endl;
  cout<<"Bonds"<<endl;
  cout<<endl;


  if(fWrite){
    outfile<<endl;
    outfile<<"Bonds"<<endl;
    outfile<<endl;
  }


  int cb = 1;
   
  for(int i = 0; i < N_bonds_tot; i++){

    if(verbose){
      if(i%2 == 0){
	cout<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
      } else {
	cout<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
      }
    }
    
    
    if(fWrite){
      if(i%2 == 0){
	outfile<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
      } else {
	outfile<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
      }
    }
    if(i%2 != 0)
      cb += 3;
    
  }
  
//#######################
//## generate angle list
//#######################

 
   cout<<endl;
   cout<<"Angles"<<endl;
   cout<<endl;
  

  if(fWrite){
    outfile<<endl;
    outfile<<"Angles"<<endl;
    outfile<<endl;
  }

  int ca = 1;
  
  //write out angle list for droplet first
  
  for(int i = 0; i < N_angles_tot; i++){

    if(verbose)
      cout<<"  "<<i+1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;
    
    if(fWrite)
      outfile<<"  "<<i+1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;
    
    ca += 3;
    
  }
  
  outfile.close();
}

void PlaceAtomAbove2ndBox(int t_step,int id){

  double dx = xb2c - xb1c;  //to shift to center of 2nd box in x direction
  double dy = yb2c - yb1c;  //to shift to center of 2nd box in y direction
  double dz = zb2hi - zb1lo; //to shift to top of 2nd box in z direction

   Atom_x1[t_step][id] = Atom_x1[t_step][id] + dx;
   Atom_y1[t_step][id] = Atom_y1[t_step][id] + dy;
   Atom_z1[t_step][id] = Atom_z1[t_step][id] + dz + dtz;  //dtz allows for separation between upper and lower box
    
}
void CutSlab(double zcent, int t_step){

  //fstream outfile("CH4onTopOfH2Ofilm_xyz.dat",ios::out);
  fstream outfile("waterslab_sapphire_xyz.dat",ios::out);
  
  N_mol1 = N_atoms1/N_atoms_mol1;
  
  //geometric box center
  double xbc1 = 0.5*(xb1lo + xb1hi);
  double ybc1 = 0.5*(yb1lo + yb1hi);
  double zbc1 = 0.5*(zb1lo + zb1hi);
  
  double dx,dy,dz,dr;
  
  int i = t_step,atid;
  int ix = 0,iy = 0,iz = 0;
  
  N_mol_sp = 0;
  
  for(int j = 0; j < N_mol1; j++){
    
    //use the center of mass of the molecule
    //    dx = Mol_xcm[i][j] - xbc;
    //    dy = Mol_ycm[i][j] - ybc;
    //    dz = Mol_zcm[i][j] - zcent;
    
    //use the oxygen atom
    atid = N_atoms_mol1*j +  1;
    dx = Atom_x1[i][atid-1] - xbc1;
    dy = Atom_y1[i][atid-1] - ybc1;
    dz = Atom_z1[i][atid-1] - zbc1;
    
    
    //    if(Atom_x1[i][atid-1] >= xSlabLo && Atom_x1[i][atid-1] <= xSlabHi &&
    //       Atom_y1[i][atid]   >= ySlabLo && Atom_y1[i][atid] <= ySlabHi &&
    //       Atom_z1[i][atid+1] >= zSlabLo && Atom_z1[i][atid+1] <= zSlabHi){
    
    //       N_mol_sp++;
    //    }    

    if(dx > xSlabLo && dx < xSlabHi &&
	dy > ySlabLo && dy < ySlabHi &&
	dz > zSlabLo && dz < zSlabHi){

      dx = Atom_x1[i][atid] - xbc1;
      dy = Atom_y1[i][atid] - ybc1;
      dz = Atom_z1[i][atid] - zbc1;

      if(dx > xSlabLo && dx < xSlabHi &&
	 dy > ySlabLo && dy < ySlabHi &&
	 dz > zSlabLo && dz < zSlabHi){

	dx = Atom_x1[i][atid+1] - xbc1;
	dy = Atom_y1[i][atid+1] - ybc1;
	dz = Atom_z1[i][atid+1] - zbc1;

	if(dx > xSlabLo && dx < xSlabHi &&
	   dy > ySlabLo && dy < ySlabHi &&
	   dz > zSlabLo && dz < zSlabHi){
		     
	  N_mol_sp++;
	}
      }
    }
  }

  int nid = 0;
  
  Atom_spid = new int[N_mol_sp*N_atoms_mol1]; 
  Atom_spmol = new int[N_mol_sp*N_atoms_mol1]; 
  Atom_sptype = new int[N_mol_sp*N_atoms_mol1]; 
  
  Atom_spx = new double[N_mol_sp*N_atoms_mol1];   
  Atom_spy = new double[N_mol_sp*N_atoms_mol1];  
  Atom_spz = new double[N_mol_sp*N_atoms_mol1];  
  
  for(int j = 0; j < N_mol1; j++){
    
    //use the center of mass of the molecule
    //    dx = Mol_xcm[i][j] - xbc;
    //    dy = Mol_ycm[i][j] - ybc;
    //    dz = Mol_zcm[i][j] - zcent;
    
    //use the oxygen atom
    atid = N_atoms_mol1*j +  1;
    dx = Atom_x1[i][atid-1] - xbc1;
    dy = Atom_y1[i][atid-1] - ybc1;
    dz = Atom_z1[i][atid-1] - zbc1;
    
    //    if(Atom_x1[i][atid-1] >= xSlabLo && Atom_x1[i][atid-1] <= xSlabHi &&
    //       Atom_y1[i][atid]   >= ySlabLo && Atom_y1[i][atid] <= ySlabHi &&
    //       Atom_z1[i][atid+1] >= zSlabLo && Atom_z1[i][atid+1] <= zSlabHi){
    if(dx > xSlabLo && dx < xSlabHi &&
       dy > ySlabLo && dy < ySlabHi &&
       dz > zSlabLo && dz < zSlabHi){
      dx = Atom_x1[i][atid] - xbc1;
      dy = Atom_y1[i][atid] - ybc1;
      dz = Atom_z1[i][atid] - zbc1;
      
      if(dx > xSlabLo && dx < xSlabHi &&
	 dy > ySlabLo && dy < ySlabHi &&
	 dz > zSlabLo && dz < zSlabHi){
	dx = Atom_x1[i][atid+1] - xbc1;
	dy = Atom_y1[i][atid+1] - ybc1;
	dz = Atom_z1[i][atid+1] - zbc1;
	
	if(dx > xSlabLo && dx < xSlabHi &&
	   dy > ySlabLo && dy < ySlabHi &&
	   dz > zSlabLo && dz < zSlabHi){
	  
	  for(int k = 0; k < 3; k++){
	    
	    id  = N_atoms_mol1*j + k + 1;
	    type = Atom_type1[i][id-1];
	    mol = Atom_mol1[i][id-1]; 
	    
	    if(fShift)
	      PlaceAtomAbove2ndBox(i,id-1);
	    
	    x = Atom_x1[i][id-1];   
	    y = Atom_y1[i][id-1];
	    z = Atom_z1[i][id-1]; 
	    
	    Atom_spid[nid] = nid;
	    Atom_spmol[nid] = j+1;
	    Atom_sptype[nid] = type;
	    Atom_spx[nid] = x;   
	    Atom_spy[nid] = y;
	    Atom_spz[nid] = z;
	    nid++; 
	    outfile<<x-xbc1<<"\t"<<y-ybc1<<"\t"<<z-zbc1<<endl;
	    if(CheckSlab)
	      cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	  }
	}
      }
    }
  }
    
  outfile.close();
}
  


void CutSphere(double zcent, int t_step){

 N_mol1 = N_atoms1/N_atoms_mol1;

 //geometric box center
 double xbc1 = 0.5*(xb1lo + xb1hi);
 double ybc1 = 0.5*(yb1lo + yb1hi);
 double zbc1 = 0.5*(zb1lo + zb1hi);
 
 
 if(Rmax == 0) Rmax = zbc1 - zb1lo - 12.0;   //avoid the edges

 // //Rmax = zbc1 - zb1lo;
//  if(zcent > zbc1)
//    Rmax = zb1hi - zcent - 2.0;
//  if(zcent < zbc1)
//    Rmax = zcent - zb1lo - 2.0;

 double dx,dy,dz,dr;

 int i = t_step,atid;
 int ix = 0,iy = 0,iz = 0;

 N_mol_sp = 0;

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
   //    dx = Mol_xcm[i][j] - xbc;
   //    dy = Mol_ycm[i][j] - ybc;
   //    dz = Mol_zcm[i][j] - zcent;
   
   //use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   dz = Atom_z1[i][atid-1] - zcent;
   
   if(fabs(dx) < Rmax && fabs(dy) < Rmax && fabs(dz) < Rmax){
     dr = sqrt(dx*dx + dy*dy + dz*dz);

     if(dr < Rmax) N_mol_sp++;
   }
 }

 int nid = 0;

 Atom_spid = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_spmol = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_sptype = new int[N_mol_sp*N_atoms_mol1]; 

 Atom_spx = new double[N_mol_sp*N_atoms_mol1];   
 Atom_spy = new double[N_mol_sp*N_atoms_mol1];  
 Atom_spz = new double[N_mol_sp*N_atoms_mol1];  

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
//    dx = Mol_xcm[i][j] - xbc;
//    dy = Mol_ycm[i][j] - ybc;
//    dz = Mol_zcm[i][j] - zcent;
   
//use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   dz = Atom_z1[i][atid-1] - zcent;
   

   if(fabs(dx) < Rmax && fabs(dy) < Rmax && fabs(dz) < Rmax){
     dr = sqrt(dx*dx + dy*dy + dz*dz);
     
     if(dr < Rmax){
       
       if(CheckDroplet)
	 cout<<"Mol: "<<j+1<<" , dx = "<<dx<<", dy = "<<dy<<" , dz = "<<dz<<", dr = "<<dr<<" , Rmax = "<<Rmax<<endl;
       
       for(int k = 0; k < 3; k++){
	 
	 id  = N_atoms_mol1*j + k + 1;
	 type = Atom_type1[i][id-1];
	 mol = Atom_mol1[i][id-1]; 
	 
	 if(fShift)
	   PlaceAtomAbove2ndBox(i,id-1);
	 
	 x = Atom_x1[i][id-1];   
	 y = Atom_y1[i][id-1];
	 z = Atom_z1[i][id-1]; 

	 Atom_spid[nid] = nid;
	 Atom_spmol[nid] = j+1;
	 Atom_sptype[nid] = type;
	 Atom_spx[nid] = x;   
	 Atom_spy[nid] = y;
	 Atom_spz[nid] = z;
	 nid++; 
	 if(CheckDroplet)
	   cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
       }
     }
   }
   
 }
}

void CutHemiSphere(double zcent, int t_step){

 N_mol1 = N_atoms1/N_atoms_mol1;

 //geometric box center
 double xbc1 = 0.5*(xb1lo + xb1hi);
 double ybc1 = 0.5*(yb1lo + yb1hi);
 double zbc1 = 0.5*(zb1lo + zb1hi);
 
 
 // if(Rmax == 0) Rmax = zbc1 - zb1lo - 12.0;   //avoid the edges

 // //Rmax = zbc1 - zb1lo;
//  if(zcent > zbc1)
//    Rmax = zb1hi - zcent - 2.0;
//  if(zcent < zbc1)
//    Rmax = zcent - zb1lo - 2.0;

 double dx,dy,dz,dr;

 int i = t_step,atid;
 int ix = 0,iy = 0,iz = 0;

 N_mol_sp = 0;

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
   //    dx = Mol_xcm[i][j] - xbc;
   //    dy = Mol_ycm[i][j] - ybc;
   //    dz = Mol_zcm[i][j] - zcent;
   
   //use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   dz = Atom_z1[i][atid-1] - zb1lo;
   
   
   if( Atom_z1[i][atid-1] < Rmax && Atom_z1[i][atid] < Rmax && Atom_z1[i][atid+1] < Rmax){
     dr = sqrt(dx*dx + dy*dy + dz*dz);
     
     if(dr < Rmax) N_mol_sp++;
   }
 }

 int nid = 0;

 Atom_spid = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_spmol = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_sptype = new int[N_mol_sp*N_atoms_mol1]; 

 Atom_spx = new double[N_mol_sp*N_atoms_mol1];   
 Atom_spy = new double[N_mol_sp*N_atoms_mol1];  
 Atom_spz = new double[N_mol_sp*N_atoms_mol1];  

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
//    dx = Mol_xcm[i][j] - xbc;
//    dy = Mol_ycm[i][j] - ybc;
//    dz = Mol_zcm[i][j] - zcent;
   
//use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   dz = Atom_z1[i][atid-1] - zb1lo;

   if( Atom_z1[i][atid-1] < Rmax && Atom_z1[i][atid] < Rmax && Atom_z1[i][atid+1] < Rmax){

     dr = sqrt(dx*dx + dy*dy + dz*dz);

     if(dr < Rmax){
       
       //if(CheckDroplet && dz > Rmax)
       if(dz < Rmax)
	 cout<<"Mol: "<<j+1<<" , dx = "<<dx<<", dy = "<<dy<<" , dz = "<<dz<<", dr = "<<dr<<" , Rmax = "<<Rmax<<endl;
       
       for(int k = 0; k < 3; k++){
	 
	 id  = atid + k;
	 type = Atom_type1[i][id-1];
	 mol = Atom_mol1[i][id-1]; 
	 
	 if(fShift)
	   PlaceAtomAbove2ndBox(i,id-1);
	 
	 x = Atom_x1[i][id-1];   
	 y = Atom_y1[i][id-1];
	 z = Atom_z1[i][id-1]; 
	 
	 Atom_spid[nid] = nid;
	 Atom_spmol[nid] = j+1;
	 Atom_sptype[nid] = type;
	 Atom_spx[nid] = x;   
	 Atom_spy[nid] = y;
	 Atom_spz[nid] = z;
	 nid++; 
	 if(nid == 12903)//CheckDroplet)
	   cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x-xbc1<<" "<<y-ybc1<<" "<<z-zb1lo<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
       }
     }
   }
   
 }
}



void CutCylinder(double zcent, int t_step){

  cout<<"In CutCylinder(double zcent, int t_step)"<<endl;
  fstream outfile(Form("spece_cylinder_%dA.dat",(int)Rmax),ios::out);

 N_mol1 = N_atoms1/N_atoms_mol1;

 //geometric box center
 double xbc1 = 0.5*(xb1lo + xb1hi);
 double ybc1 = 0.5*(yb1lo + yb1hi);
 double zbc1 = 0.5*(zb1lo + zb1hi);
 
 if(zcent > zb1hi || zcent < zb1lo) {
   cout<<"error: zcent outside of box"<<endl;
   return;
 }

 if(zcent > zb1hi || zcent < zb1lo) {
   cout<<"error: zcent outside of box"<<endl;
   return;
 }
//  if(zcent > zbc1 && Rmax > zb1hi-zcent)
//    Rmax = fabs(zb1hi-zcent-2.0);
//  if(zcent < zbc1 && Rmax > zcent-zb1lo)
//    Rmax = fabs(zcent-zb1lo-2.0);
 // if(Rmax == 0) Rmax = zbc1 - zb1lo - 12.0;   //avoid the edges

 // //Rmax = zbc1 - zb1lo;
//  if(zcent > zbc1)
//    Rmax = zb1hi - zcent - 2.0;
//  if(zcent < zbc1)
//    Rmax = zcent - zb1lo - 2.0;

 double dx,dy,dz,dr;

 int i = t_step,atid;
 int ix = 0,iy = 0,iz = 0;

 N_mol_sp = 0;

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
   //    dx = Mol_xcm[i][j] - xbc;
   //    dy = Mol_ycm[i][j] - ybc;
   //    dz = Mol_zcm[i][j] - zcent;
   
   //use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   //dz = Atom_z1[i][atid-1] - zcent;
   dz = Atom_z1[i][atid-1] - zbc1;
   
//    if(fabs(dy) < 0.5*yCylinderLength){
//      dr = sqrt(dx*dx + dz*dz);
//      if(dr < Rmax) N_mol_sp++;
//    }


   if(Atom_y1[i][atid-1] >= yCylinderLo && Atom_y1[i][atid-1] <= yCylinderHi &&
      Atom_y1[i][atid] >= yCylinderLo && Atom_y1[i][atid] <= yCylinderHi &&
      Atom_y1[i][atid+1] >= yCylinderLo && Atom_y1[i][atid+1] <= yCylinderHi){
     
     dr = sqrt(dx*dx + dz*dz);
     if(dr < Rmax) N_mol_sp++;
   }
   
   
 }

 int nid = 0;

 Atom_spid = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_spmol = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_sptype = new int[N_mol_sp*N_atoms_mol1]; 

 Atom_spx = new double[N_mol_sp*N_atoms_mol1];   
 Atom_spy = new double[N_mol_sp*N_atoms_mol1];  
 Atom_spz = new double[N_mol_sp*N_atoms_mol1];  

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
//    dx = Mol_xcm[i][j] - xbc;
//    dy = Mol_ycm[i][j] - ybc;
//    dz = Mol_zcm[i][j] - zcent;
   
//use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   //dz = Atom_z1[i][atid-1] - zcent;
   dz = Atom_z1[i][atid-1] - zbc1;


   if(Atom_y1[i][atid-1] >= yCylinderLo && Atom_y1[i][atid-1] <= yCylinderHi &&
      Atom_y1[i][atid] >= yCylinderLo && Atom_y1[i][atid] <= yCylinderHi &&
      Atom_y1[i][atid+1] >= yCylinderLo && Atom_y1[i][atid+1] <= yCylinderHi){
     
     dr = sqrt(dx*dx + dz*dz);
  
     if(dr < Rmax){
       
       //if(CheckDroplet && dz > Rmax)
//        if(dz < Rmax)
// 	 cout<<"Mol: "<<j+1<<" , dx = "<<dx<<", dy = "<<dy<<" , dz = "<<dz<<", dr = "<<dr<<" , Rmax = "<<Rmax<<endl;
       
       for(int k = 0; k < 3; k++){
	 
	 id  = atid + k;
	 type = Atom_type1[i][id-1];
	 mol = Atom_mol1[i][id-1]; 
	 
	 if(fShift)
	   PlaceAtomAbove2ndBox(i,id-1);
	 
	 x = Atom_x1[i][id-1];   
	 y = Atom_y1[i][id-1];
	 z = Atom_z1[i][id-1]; 
	 
	 Atom_spid[nid] = nid;
	 Atom_spmol[nid] = j+1;
	 Atom_sptype[nid] = type;
	 Atom_spx[nid] = x;   
	 Atom_spy[nid] = y;
	 Atom_spz[nid] = z;

	 x -= xbc1;
	 z -= zbc1;

	 nid++; 
// 	 if(nid == 12903)//CheckDroplet)
// 	   cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x-xbc1<<" "<<y-ybc1<<" "<<z-zb1lo<<" "<<ix<<" "<<iy<<" "<<iz<<endl;

	 if(Atom_y1[i][atid-1] >= yCylinderLo && Atom_y1[i][atid-1] <= yCylinderHi){
	   //cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x-xbc1<<" "<<y-ybc1<<" "<<z-zb1lo<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	   //cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	   //cout<<x<<"\t"<<y<<"\t"<<z<<endl;
	   outfile<<x<<"\t"<<y<<"\t"<<z<<endl;
	 }
       }
     }
   }
   
 }
}

void CutHemiCylinder(double zcent, int t_step){

  cout<<"In CutHemiCylinder(double zcent, int t_step)"<<endl;
  fstream outfile(Form("spece_hemicylinder_%dA.dat",(int)Rmax),ios::out);


 N_mol1 = N_atoms1/N_atoms_mol1;

 //geometric box center
 double xbc1 = 0.5*(xb1lo + xb1hi);
 double ybc1 = 0.5*(yb1lo + yb1hi);
 double zbc1 = 0.5*(zb1lo + zb1hi);
 
 cout<<"xb1c = "<<xbc1<<" , yb1c = "<<ybc1<<" , zb1c = "<<zbc1<<endl;

 if(zcent > zb1hi || zcent < zb1lo) {
   cout<<"error: zcent outside of box"<<endl;
   return;
 }

 zcent = zb1lo + 10;

//  if(zcent > zbc1 && Rmax > zb1hi-zcent)
//    Rmax = fabs(zb1hi-zcent-2.0);
//  if(zcent < zbc1 && Rmax > zcent-zb1lo)
//    Rmax = fabs(zcent-zb1lo-2.0);
 // if(Rmax == 0) Rmax = zbc1 - zb1lo - 12.0;   //avoid the edges

 // //Rmax = zbc1 - zb1lo;
//  if(zcent > zbc1)
//    Rmax = zb1hi - zcent - 2.0;
//  if(zcent < zbc1)
//    Rmax = zcent - zb1lo - 2.0;

 double dx,dy,dz,dr;
 double spx,spy,xH1,yH1,xH2,yH2;

 int i = t_step,atid;
 int ix = 0,iy = 0,iz = 0;

 N_mol_sp = 0;

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
   //    dx = Mol_xcm[i][j] - xbc;
   //    dy = Mol_ycm[i][j] - ybc;
   //    dz = Mol_zcm[i][j] - zcent;
   
   //use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   spx = Atom_x1[i][atid-1];
   spy = Atom_y1[i][atid-1];
   xH1 =  Atom_x1[i][atid]; 
   yH1 =  Atom_y1[i][atid]; 
   xH2 =  Atom_x1[i][atid+1]; 
   yH2 =  Atom_y1[i][atid+1];     
 
   if(TMath::Abs(xH1 - spx) > 0.5*Lx){
     if(xH1 < spx)
       Atom_x1[i][atid] = xH1 + Lx;
     else 
       Atom_x1[i][atid] = xH1 - Lx;
   }
   if(TMath::Abs(xH2 - spx) > 0.5*Lx){
     if(xH2 < spx)
       Atom_x1[i][atid+1] = xH2 + Lx;
     else 
       Atom_x1[i][atid+1] = xH2 - Lx; 
   } 
   if(TMath::Abs(yH1 - spy) > 0.5*Ly){
     if(yH1 < spy)
       Atom_y1[i][atid] = yH1 + Ly;
     else 
       Atom_y1[i][atid] = yH1 - Ly;
   }
   if(TMath::Abs(yH2 - spy) > 0.5*Ly){
     if(yH2 < spy)
       Atom_y1[i][atid+1] = yH2 + Ly;
     else 
       Atom_y1[i][atid+1] = yH2 - Ly;
   }
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   dz = Atom_z1[i][atid-1] - zcent;
 
//    if(fabs(dy) < 0.5*yCylinderLength && Atom_z1[i][atid-1] > zcent){
//      dr = sqrt(dx*dx + dz*dz);
//      if(dr < Rmax) N_mol_sp++;
//    }

   if(Atom_z1[i][atid-1] >= zcent){
     if(Atom_y1[i][atid-1] >= yCylinderLo && Atom_y1[i][atid-1] <= yCylinderHi &&
	Atom_y1[i][atid] >= yCylinderLo && Atom_y1[i][atid] <= yCylinderHi &&
	Atom_y1[i][atid+1] >= yCylinderLo && Atom_y1[i][atid+1] <= yCylinderHi){
       dr = sqrt(dx*dx + dz*dz);
       if(dr < Rmax) N_mol_sp++;
     }
   }
 }

 cout<<"N_mol_sp = "<<N_mol_sp<<endl;
 int nid = 0;

 Atom_spid = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_spmol = new int[N_mol_sp*N_atoms_mol1]; 
 Atom_sptype = new int[N_mol_sp*N_atoms_mol1]; 

 Atom_spx = new double[N_mol_sp*N_atoms_mol1];   
 Atom_spy = new double[N_mol_sp*N_atoms_mol1];  
 Atom_spz = new double[N_mol_sp*N_atoms_mol1];  

 for(int j = 0; j < N_mol1; j++){
   
   //use the center of mass of the molecule
//    dx = Mol_xcm[i][j] - xbc;
//    dy = Mol_ycm[i][j] - ybc;
//    dz = Mol_zcm[i][j] - zcent;
   
//use the oxygen atom
   atid = N_atoms_mol1*j +  1;
   dx = Atom_x1[i][atid-1] - xbc1;
   dy = Atom_y1[i][atid-1] - ybc1;
   dz = Atom_z1[i][atid-1] - zcent;


   if(Atom_z1[i][atid-1] >= zcent){
     if(Atom_y1[i][atid-1] >= yCylinderLo && Atom_y1[i][atid-1] <= yCylinderHi &&
	Atom_y1[i][atid] >= yCylinderLo && Atom_y1[i][atid] <= yCylinderHi &&
	Atom_y1[i][atid+1] >= yCylinderLo && Atom_y1[i][atid+1] <= yCylinderHi){
       
       dr = sqrt(dx*dx + dz*dz);
       
       if(dr < Rmax){
	 
	 //if(CheckDroplet && dz > Rmax)
	 //if(dz < Rmax)
	 //cout<<"Mol: "<<j+1<<" , dx = "<<dx<<", dy = "<<dy<<" , dz = "<<dz<<", dr = "<<dr<<" , Rmax = "<<Rmax<<endl;
	 
	 for(int k = 0; k < 3; k++){
	   
	   id  = atid + k;
	   type = Atom_type1[i][id-1];
	   mol = Atom_mol1[i][id-1]; 
	   
	   if(fShift)
	     PlaceAtomAbove2ndBox(i,id-1);
	   
	   x = Atom_x1[i][id-1];   
	   y = Atom_y1[i][id-1];
	   z = Atom_z1[i][id-1]; 
	   
	   Atom_spid[nid] = nid;
	   Atom_spmol[nid] = j+1;
	   Atom_sptype[nid] = type;
	   Atom_spx[nid] = x;   
	   Atom_spy[nid] = y;
	   Atom_spz[nid] = z;
	   nid++; 
	   //if(nid == 12903)//CheckDroplet)
	   //if(Atom_y1[i][atid-1] <= yCylinderLo || Atom_y1[i][atid-1] >= yCylinderHi){
	   if(Atom_y1[i][atid-1] >= yCylinderLo && Atom_y1[i][atid-1] <= yCylinderHi){
	     //cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x-xbc1<<" "<<y-ybc1<<" "<<z-zb1lo<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	     //cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	     //cout<<x<<"\t"<<y<<"\t"<<z<<endl;
	     outfile<<x<<"\t"<<y<<"\t"<<z<<endl;
	   }
	 }
       }
     }
     
   }
 }


 outfile.close();

 cout<<"Done in CutHemiCylinder(double zcent, int t_step)"<<endl;
}


void Initialize(){

  M = Mass[0] + 2.*Mass[1];  //mass/mole of a water (in grams)

  Time_step1 = new int[n_time_steps1];
  Atom_id1   = new int*[n_time_steps1];
  Atom_type1 = new int*[n_time_steps1];
  Atom_mol1  = new int*[n_time_steps1];

  Atom_x1    = new double*[n_time_steps1];
  Atom_y1    = new double*[n_time_steps1];
  Atom_z1    = new double*[n_time_steps1];

  Vol1 = new double[n_time_steps1];
  rho1 = new double[n_time_steps1];
  N_den1 = new double[n_time_steps1];

  Mol_xcm1    = new double*[n_time_steps1];
  Mol_ycm1    = new double*[n_time_steps1];
  Mol_zcm1    = new double*[n_time_steps1];

  Time_step2 = new int[n_time_steps2];
  Atom_id2   = new int*[n_time_steps2];
  Atom_type2 = new int*[n_time_steps2];
  Atom_mol2  = new int*[n_time_steps2];

  Atom_x2    = new double*[n_time_steps2];
  Atom_y2    = new double*[n_time_steps2];
  Atom_z2    = new double*[n_time_steps2];

  Vol2 = new double[n_time_steps2];
  rho2 = new double[n_time_steps2];
  N_den2 = new double[n_time_steps2];

  Mol_xcm2    = new double*[n_time_steps2];
  Mol_ycm2    = new double*[n_time_steps2];
  Mol_zcm2    = new double*[n_time_steps2];
}


void ReadFile1(){


  //### declare variables for reading

// ITEM: TIMESTEP
// 0
// ITEM: NUMBER OF ATOMS
// 3993
// ITEM: BOX BOUNDS pp pp pp
// 0 33.723
// 0 33.723
// 0 33.723
// ITEM: ATOMS id type mol x y z vx vy vz 
// 1 1 1 1.53286 1.53286 1.53286 -0.000836477 -0.000681187 0.00172607 

  //### open file

  char line[250];

  //ifstream infile("dumpwaternpt626262_2.water",ios::in);
  ifstream infile("dump_lastTimeStep.dat",ios::in);
  //ifstream infile("waterdata120120120.spce",ios::in);
  //ifstream infile("atom_npt_c.dat",ios::in);   //for PMMA for Mesfin

  if (infile.is_open())
    {

      for(int i = 0; i < n_time_steps1; i++){
	//while ( infile.good() )
	//read first 9 lines
	
	infile.getline(line,250);
	cout << line << endl;
	infile>>time_step;
	cout<<time_step<<endl;
	infile.ignore();
	infile.getline(line,250);
	cout << line << endl;
	infile>>N_atoms1;
	cout<<N_atoms1<<endl;
	infile.ignore();
	infile.getline(line,250);
	cout << line << endl;	 
	infile>>xlo>>xhi>>ylo>>yhi>>zlo>>zhi;
	cout<<xlo<<" "<<xhi<<endl;
	cout<<ylo<<" "<<yhi<<endl;
	cout<<zlo<<" "<<zhi<<endl;
	infile.ignore();
	infile.getline(line,250);
	cout << line << endl;
     
	xb1lo = xlo;  xb1hi = xhi;
	yb1lo = ylo;  yb1hi = yhi;
	zb1lo = zlo;  zb1hi = zhi;

// 	//calculate box center
// 	if(i == 0)
// 	  CalcBoxCenter();

	//compute volume;
	Lx = xhi-xlo;
	Ly = yhi-ylo;
	Lz = zhi-zlo;	

	Vol1[i] = Lx*Ly*Lz;

	N_mol1 = N_atoms1/N_atoms_mol1;
	//compute density;
	rho1[i] = ((N_mol1*M/NA)/Vol1[i])/rho_exact; 
	N_den1[i] = N_mol1/Vol1[i]; 


	Time_step1[i] = time_step;
	Atom_id1[i]   = new int[N_atoms1];
	Atom_type1[i] = new int[N_atoms1];
	Atom_mol1[i]  = new int[N_atoms1];

	Atom_x1[i]    = new double[N_atoms1];
	Atom_y1[i]    = new double[N_atoms1];
	Atom_z1[i]    = new double[N_atoms1];
	
	Mol_xcm1[i]   = new double[N_mol1];
	Mol_ycm1[i]   = new double[N_mol1];
	Mol_zcm1[i]   = new double[N_mol1];

  //for center of mass calculation, zero the array element first
	 for(int ml = 0; ml < N_mol1; ml++){
	    Mol_xcm1[i][ml] = 0.;
	    Mol_ycm1[i][ml] = 0.;
	    Mol_zcm1[i][ml] = 0.;
	 }
	//read atom information

	//N_atoms = 10;
	
	for(int k = 0; k < N_atoms1; k++){


	  infile>>id>>type>>mol>>xs>>ys>>zs>>ix>>iy>>iz;

	  Atom_id1[i][id-1]   = id;
	  Atom_type1[i][id-1] = type;
	  Atom_mol1[i][id-1]  = mol; 

	  //actual position
	  xu = (xs + ix)*Lx;   
	  yu = (ys + iy)*Ly;
	  zu = (zs + iz)*Lz; 
	  //position wrapped into box 

	    //position wrapped Int_to box 
// 	  x = xlo + xs*Lx;   
// 	  y = ylo + ys*Ly;
// 	  z = zlo + zs*Lz;

	  x = xs*Lx;   
	  y = ys*Ly;
	  z = zs*Lz;
 
	  if(UseActualPos){
	    Atom_x1[i][id-1] = xu;
	    Atom_y1[i][id-1] = yu;
	    Atom_z1[i][id-1] = zu;
	  } else {
	    Atom_x1[i][id-1] = xu;
	    Atom_y1[i][id-1] = yu;
	    Atom_z1[i][id-1] = zu;
// 	    Atom_x1[i][id-1] = x;   
// 	    Atom_y1[i][id-1] = y;
// 	    Atom_z1[i][id-1] = z;
	    
	  } 
 

	  //sum m*x etc for numerator of CM expression and divide outside this loop by the total mass of the molecule
	  cout << "Atom type = " << type << endl;
	  Mol_xcm1[i][mol-1] += Mass[type - 1]*Atom_x1[i][id-1];  
	  Mol_ycm1[i][mol-1] += Mass[type - 1]*Atom_y1[i][id-1];  
	  Mol_zcm1[i][mol-1] += Mass[type - 1]*Atom_z1[i][id-1];  

// 	  if(UseActualPos)
// 	    cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (xu,yu,zu,ix,iy,iz): "<<xu<<" "<<yu<<" "<<zu<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
// 	  else
// 	    cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z,ix,iy,iz): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
	  
	}
	//divide by the total mass of the molecule
	for(int ml = 0; ml < N_mol1; ml++){
	  Mol_xcm1[i][ml] /= M;
	  Mol_ycm1[i][ml] /= M;
	  Mol_zcm1[i][ml] /= M;
	}

	//infile.getline(line,250);
	infile.ignore();
	infile.getline(line,250);
      }
      infile.close();
    }


  if(Check){
    for(int i = n_time_steps1 - 1; i < n_time_steps1; i++){
      for(int k = 0; k < N_atoms1; k++){
	
	id = Atom_id1[i][k];
	type = Atom_type1[i][k];
	mol = Atom_mol1[i][k];
	x = Atom_x1[i][k];   
	y = Atom_y1[i][k];
	z = Atom_z1[i][k];

	cout<<"Time_step; "<<i*dump_steps1<<" , atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<endl; 
      }
    }
  }

  N_atoms1 = N_mol1*N_atoms_mol1;
  N_bonds1 = N_mol1*Nbonds_mol1;
  N_angles1 = N_mol1*Nangles_mol1;
  N_dihedrals1 = 0;
  N_impropers1 = 0;
}

void ReadFile2(){


  //### declare variables for reading

// ITEM: TIMESTEP
// 0
// ITEM: NUMBER OF ATOMS
// 3993
// ITEM: BOX BOUNDS pp pp pp
// 0 33.723
// 0 33.723
// 0 33.723
// ITEM: ATOMS id type mol x y z vx vy vz 
// 1 1 1 1.53286 1.53286 1.53286 -0.000836477 -0.000681187 0.00172607 


  //### open file

  char line[250];

  //ifstream infile("dumpwaternpt626262_2.water",ios::in);
  
  ifstream infile("dump_lastTimeStep.dat",ios::in);
  //ifstream infile("atom_npt_c.dat",ios::in);   //for PMMA for Mesfin

  if (infile.is_open())
    {

      for(int i = 0; i < n_time_steps2; i++){
	//while ( infile.good() )
	//read first 9 lines
	
	infile.getline(line,250);
	cout << line << endl;
	infile>>time_step;
	cout<<time_step<<endl;
	infile.ignore();
	infile.getline(line,250);
	cout << line << endl;
	infile>>N_atoms2;
	cout<<N_atoms2<<endl;
	infile.ignore();
	infile.getline(line,250);
	cout << line << endl;	 
	infile>>xlo>>xhi>>ylo>>yhi>>zlo>>zhi;
	cout<<xlo<<" "<<xhi<<endl;
	cout<<ylo<<" "<<yhi<<endl;
	cout<<zlo<<" "<<zhi<<endl;
	infile.ignore();
	infile.getline(line,250);
	cout << line << endl;
     
	xb2lo = xlo;  xb2hi = xhi;
	yb2lo = ylo;  yb2hi = yhi;
	zb2lo = zlo;  zb2hi = zhi;

// 	//calculate box center
// 	if(i == 0)
// 	  CalcBoxCenter();

	//compute volume;
	Lx = xhi-xlo;
	Ly = yhi-ylo;
	Lz = zhi-zlo;	

	Vol2[i] = Lx*Ly*Lz;

	N_mol2 = N_atoms2/N_atoms_mol2;
	//compute density;
	rho2[i] = ((N_mol2*M/NA)/Vol2[i])/rho_exact; 
	N_den2[i] = N_mol2/Vol2[i]; 


	Time_step2[i] = time_step;
	Atom_id2[i]   = new int[N_atoms2];
	Atom_type2[i] = new int[N_atoms2];
	Atom_mol2[i]  = new int[N_atoms2];

	Atom_x2[i]    = new double[N_atoms2];
	Atom_y2[i]    = new double[N_atoms2];
	Atom_z2[i]    = new double[N_atoms2];

	Mol_xcm2[i]   = new double[N_mol2];
	Mol_ycm2[i]   = new double[N_mol2];
	Mol_zcm2[i]   = new double[N_mol2];

  //for center of mass calculation, zero the array element first
	 for(int ml = 0; ml < N_mol2; ml++){
	    Mol_xcm2[i][ml] = 0.;
	    Mol_ycm2[i][ml] = 0.;
	    Mol_zcm2[i][ml] = 0.;
	 }
	//read atom information

	//N_atoms = 10;
	
	for(int k = 0; k < N_atoms2; k++){

	  
	  infile>>id>>type>>mol>>xs>>ys>>zs>>ix>>iy>>iz;
 

	  Atom_id2[i][id-1]   = id;
	  Atom_type2[i][id-1] = type;
	  Atom_mol2[i][id-1]  = mol; 

	  //actual position
	  xu = (xs + ix)*Lx;   
	  yu = (ys + iy)*Ly;
	  zu = (zs + iz)*Lz; 
	  //position wrapped into box 
	  x = xs*Lx;   
	  y = ys*Ly;
	  z = zs*Lz;
 
	  if(UseActualPos){
	    Atom_x2[i][id-1] = xu;
	    Atom_y2[i][id-1] = yu;
	    Atom_z2[i][id-1] = zu;
	  } else {
	    Atom_x2[i][id-1] = x;   
	    Atom_y2[i][id-1] = y;
	    Atom_z2[i][id-1] = z; 
	  } 
 
	
	  //sum m*x etc for numerator of CM expression and divide outside this loop by the total mass of the molecule
	  Mol_xcm2[i][mol-1] += Mass[type - 1]*Atom_x2[i][id-1];  
	  Mol_ycm2[i][mol-1] += Mass[type - 1]*Atom_y2[i][id-1];  
	  Mol_zcm2[i][mol-1] += Mass[type - 1]*Atom_z2[i][id-1];  

// 	  if(UseActualPos)
// 	    cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (xu,yu,zu,ix,iy,iz): "<<xu<<" "<<yu<<" "<<zu<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
// 	  else
// 	    cout<<"atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z,ix,iy,iz): "<<x<<" "<<y<<" "<<z<<" "<<ix<<" "<<iy<<" "<<iz<<endl;

	}
	//divide by the total mass of the molecule
	for(int ml = 0; ml < N_mol2; ml++){
	  Mol_xcm2[i][ml] /= M;
	  Mol_ycm2[i][ml] /= M;
	  Mol_zcm2[i][ml] /= M;
	}

	//infile.getline(line,250);
	infile.ignore();
	infile.getline(line,250);
      }
      infile.close();
    }


  if(Check){
    for(int i = n_time_steps2 - 2; i < n_time_steps2; i++){
      for(int k = 0; k < N_atoms2; k++){
	
	id = Atom_id2[i][k];
	type = Atom_type2[i][k];
	mol = Atom_mol2[i][k];
	x = Atom_x2[i][k];   
	y = Atom_y2[i][k];
	z = Atom_z2[i][k];
	
	cout<<"Time_step; "<<i*dump_steps2<<" , atom id: "<<id<<" , atom type: "<<type<<" , Mol id:"<<mol<<" , (x,y,z): "<<x<<" "<<y<<" "<<z<<endl; 
      }
    }
  }


  N_atoms2 = N_mol2*N_atoms_mol2;
  N_bonds2 = N_mol2*Nbonds_mol2;
  N_angles2 = N_mol2*Nangles_mol2;
  N_dihedrals2 = 0;
  N_impropers2 = 0;

}

void CombineList(int t_step){

  fstream outfile("system_full_626262.spce",ios::out);

  int N_atoms_sp = N_mol_sp*N_atoms_mol1;
  int N_bonds_sp = N_mol_sp*Nbonds_mol1;
  int N_angles_sp = N_mol_sp*Nangles_mol1;
  int N_dihedrals_sp = 0;
  int N_impropers_sp = 0;

  int N_atoms_tot = N_atoms_sp + N_atoms2;
  int N_bonds_tot = N_bonds_sp + N_bonds2;
  int N_angles_tot = N_angles_sp + N_angles2;
  int N_dihedrals_tot = N_dihedrals_sp +  N_dihedrals2;
  int N_impropers_tot = N_impropers_sp + N_impropers2;

  cout<<"LAMMPS Atom File"<<endl;
  cout<<endl;
  cout<<"          "<<N_atoms_tot<<" atoms"<<endl;
  cout<<"          "<<N_bonds_tot<<" bonds"<<endl;
  cout<<"          "<<N_angles_tot<<" angles"<<endl;
  cout<<"          "<<N_dihedrals_tot<<" dihedrals"<<endl;
  cout<<"          "<<N_impropers_tot<<" impropers"<<endl;
  cout<<endl;
  cout<<"          "<<N_atom_types1<<" atom types"<<endl;
  cout<<"          "<<N_bond_types1<<" bond types"<<endl;
  cout<<"          "<<N_angle_types1<<" angle types"<<endl;
  cout<<endl;
  
  cout<<"     "<<xb2lo<<"    "<<xb2hi<<"  xlo xhi"<<endl;
  cout<<"     "<<yb2lo<<"    "<<yb2hi<<"  ylo yhi"<<endl;
  cout<<"     "<<zb2lo<<"    "<<zb2hi + (zb1hi - zb1lo) <<"  zlo zhi"<<endl;
  
  cout<<endl;  
  cout<<"Masses"<<endl;
  cout<<endl;
  cout<<"\t"<<1<<"\t"<<Mass[0]<<endl;
  cout<<"\t"<<2<<"\t"<<Mass[1]<<endl;
  

  if(fWrite){
   
    outfile<<"LAMMPS Atom File"<<endl;
    outfile<<endl;
    outfile<<"          "<<N_atoms_tot<<" atoms"<<endl;
    outfile<<"          "<<N_bonds_tot<<" bonds"<<endl;
    outfile<<"          "<<N_angles_tot<<" angles"<<endl;
    outfile<<"          "<<N_dihedrals_tot<<" dihedrals"<<endl;
    outfile<<"          "<<N_impropers_tot<<" impropers"<<endl;
    outfile<<endl;
    outfile<<"          "<<N_atom_types1<<" atom types"<<endl;
    outfile<<"          "<<N_bond_types1<<" bond types"<<endl;
    outfile<<"          "<<N_angle_types1<<" angle types"<<endl;
    outfile<<endl;
    
    outfile<<"     "<<xb2lo<<"    "<<xb2hi<<"  xlo xhi"<<endl;
    outfile<<"     "<<yb2lo<<"    "<<yb2hi<<"  ylo yhi"<<endl;
    outfile<<"     "<<zb2lo<<"    "<<zb2hi + (zb1hi - zb1lo) <<"  zlo zhi"<<endl;
  
    outfile<<endl;  
    outfile<<"Masses"<<endl;
    outfile<<endl;
    outfile<<"\t"<<1<<"\t"<<Mass[0]<<endl;
    outfile<<"\t"<<2<<"\t"<<Mass[1]<<endl;
    
 
  }
  
  //##############################
  //## generate atom list
  //##############################
	  
    cout<<endl;
    cout<<"Atoms"<<endl;
    cout<<endl;

  if(fWrite){
    outfile<<endl;
    outfile<<"Atoms"<<endl;
    outfile<<endl;
  }

  int nid,mid = 0,type;
  double ch,x,y,z;


  //write out droplet first
  for(int m = 0; m < N_atoms_sp; m++){
    
    //nid  = Atom_spid[m];
    //mid  = Atom_spmol[m];
    nid = m+1;
    if(m%N_atoms_mol1 == 0) mid++;
    type = Atom_sptype[m];
    x  = Atom_spx[m];   
    y  = Atom_spy[m];
    z  = Atom_spz[m];

    if(m%3 == 0)
      ch = charge[0];
    else
      ch = charge[1];

    cout<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	<<x<<"  "<<y<<"  "<<z<<"  0  0  0"<<endl;	
    outfile<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	<<x<<"  "<<y<<"  "<<z<<"  0  0  0"<<endl;	
  
//       1    1  1 -0.8472   12.12456   28.09298   22.27452  0  1  0
//       2    1  2  0.4236   12.53683   28.75606   22.89928  0  1  0
//       3    1  2  0.4236   11.49482   28.56390   21.65678  0  1  0
  }


  for(int m = 0; m < N_atoms2; m++){
    
    //nid  = Atom_spid[m];
    //mid  = Atom_spmol[m];
    nid = N_atoms_sp + m+1;
    if(m%N_atoms_mol1 == 0) mid++;
    type = Atom_type2[t_step][m];
    x  = Atom_x2[t_step][m];   
    y  = Atom_y2[t_step][m];
    z  = Atom_z2[t_step][m];

    if(m%3 == 0)
      ch = charge[0];
    else
      ch = charge[1];

    cout<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	<<x<<"  "<<y<<"  "<<z<<"  0  0  0"<<endl;	
    outfile<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	<<x<<"  "<<y<<"  "<<z<<"  0  0  0"<<endl;	
  
//       1    1  1 -0.8472   12.12456   28.09298   22.27452  0  1  0
//       2    1  2  0.4236   12.53683   28.75606   22.89928  0  1  0
//       3    1  2  0.4236   11.49482   28.56390   21.65678  0  1  0
  }


// //########################
// //## generate bond list
// //########################
  

    cout<<endl;
    cout<<"Bonds"<<endl;
    cout<<endl;


  if(fWrite){
    outfile<<endl;
    outfile<<"Bonds"<<endl;
    outfile<<endl;
  }

   double cb = 1;
   //write out bond list for droplet first
   for(int i = 0; i < N_bonds_sp; i++){

  
     if(i%2 == 0){
       cout<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
     } else {
	cout<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
     }
   

    if(fWrite){
      if(i%2 == 0){
	outfile<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
      } else {
	outfile<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
      }
    }
    if(i%2 != 0)
      cb += 3;
    
  }

   int nb;
   //write out bond list for substrate
  for(int i = 0; i < N_bonds2; i++){

    nb = N_bonds_sp + i;
    if(nb%2 == 0){
      cout<<"  "<<nb + 1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
    } else {
      cout<<"  "<<nb + 1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
    }
    
    
    if(fWrite){
      if(nb%2 == 0){
	outfile<<"  "<<nb + 1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
      } else {
	outfile<<"  "<<nb + 1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
      }
    }
    if(nb%2 != 0)
      cb += 3;
    
  }


   //write out bond list for substrate


//       1   1      1      2
//       2   1      1      3
//       3   1      4      5
//       4   1      4      6

//#######################
//## generate angle list
//#######################

 
   cout<<endl;
   cout<<"Angles"<<endl;
   cout<<endl;
  

  if(fWrite){
    outfile<<endl;
    outfile<<"Angles"<<endl;
    outfile<<endl;
  }

  double ca = 1;

  //write out angle list for droplet first

  for(int i = 0; i < N_angles_sp; i++){

    cout<<"  "<<i+1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;

    if(fWrite)
      outfile<<"  "<<i+1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;
    
    ca += 3;
    
  }

  //write out angle list for substrate

  int na;
  for(int i = 0; i < N_angles2; i++){

    na = N_angles_sp + i;
    cout<<"  "<<na + 1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;

    if(fWrite)
      outfile<<"  "<<na + 1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;
    
    ca += 3;
    
  }


//       1   1      2      1      3
//       2   1      5      4      6
//       3   1      8      7      9
//       4   1     11     10     12
//       5   1     14     13     15
//       6   1     17     16     18

   outfile.close();

}

void MakeSphereList(TString name){


  cout<<"########################################"<<endl;
  cout<<"Writing data to file ... "<<endl;
  cout<<"########################################"<<endl;

  fstream outfile;

  //outfile.open(name.c_str(),ios::out);
  outfile.open(name.Data(),ios::out);

  //geometric box center
  double xbc1 = 0.5*(xb1lo + xb1hi);
  double ybc1 = 0.5*(yb1lo + yb1hi);
  double zbc1 = 0.5*(zb1lo + zb1hi);
  xb1c = xbc1;
  yb1c = ybc1;
  zb1c = zbc1;

//   xb1lo = xb1lo - xb1c;
//   xb1hi = xb1hi - xb1c;
//   yb1lo = yb1lo - yb1c;
//   yb1hi = yb1hi - yb1c;
//   zb1lo = zb1lo - zb1c;
//   zb1hi = zb1hi - zb1c;

  xb1lo = xSlabLo;
  xb1hi = xSlabHi;
  yb1lo = ySlabLo;
  yb1hi = ySlabHi;
  zb1lo = zSlabLo;
  zb1hi = zSlabHi;

  cout<<"center of duplicated box"<<endl;
  cout<<"xb1c = "<<xb1c<<" , yb1c = "<<yb1c<<" , zb1c = "<<zb1c<<endl;

  cout<<"N_mol_sp: "<<N_mol_sp<<endl;
  N_atoms_sp = N_mol_sp*N_atoms_mol1;
  N_bonds_sp = N_mol_sp*Nbonds_mol1;
  N_angles_sp = N_mol_sp*Nangles_mol1;
  N_dihedrals_sp = 0;
  N_impropers_sp = 0;

  if(CheckDroplet){
    cout<<"LAMMPS Atom File"<<endl;
    cout<<endl;
    cout<<"          "<<N_atoms_sp<<" atoms"<<endl;
    cout<<"          "<<N_bonds_sp<<" bonds"<<endl;
    cout<<"          "<<N_angles_sp<<" angles"<<endl;
    cout<<"          "<<N_dihedrals_sp<<" dihedrals"<<endl;
    cout<<"          "<<N_impropers_sp<<" impropers"<<endl;
    cout<<endl;
    cout<<"          "<<N_atom_types1<<" atom types"<<endl;
    cout<<"          "<<N_bond_types1<<" bond types"<<endl;
    cout<<"          "<<N_angle_types1<<" angle types"<<endl;
    cout<<endl;
    
    if(fCombine){
      cout<<"     "<<xb2lo<<"    "<<xb2hi<<"  xlo xhi"<<endl;
      cout<<"     "<<yb2lo<<"    "<<yb2hi<<"  ylo yhi"<<endl;
      
      if(fShift)
	cout<<"     "<<zb2lo<<"    "<<zb2hi + (zb1hi - zb1lo) <<"  zlo zhi"<<endl;
      else
	cout<<"     "<<zb2lo<<"    "<<zb2hi<<"  zlo zhi"<<endl;

    } else {
      cout<<"     "<<xb1lo-xb1c<<"    "<<xb1hi-xb1c<<"  xlo xhi"<<endl;
      cout<<"     "<<yb1lo-yb1c<<"    "<<yb1hi-yb1c<<"  ylo yhi"<<endl;
      
      if(fShift)
	cout<<"     "<<zb1lo-zb1c<<"    "<<zb1hi-zb1c<<"  zlo zhi"<<endl;
      else
	cout<<"     "<<zb1lo-zb1c<<"    "<<zb1hi-zb1c<<"  zlo zhi"<<endl;
   
    }

    cout<<endl;  
    cout<<"Masses"<<endl;
    cout<<endl;
    cout<<"\t"<<1<<"\t"<<Mass[0]<<endl;
    cout<<"\t"<<2<<"\t"<<Mass[1]<<endl;
  }

  if(fWrite){
   
    outfile<<"LAMMPS Atom File"<<endl;
    outfile<<endl;
    outfile<<"          "<<N_atoms_sp<<" atoms"<<endl;
    outfile<<"          "<<N_bonds_sp<<" bonds"<<endl;
    outfile<<"          "<<N_angles_sp<<" angles"<<endl;
    outfile<<"          "<<N_dihedrals_sp<<" dihedrals"<<endl;
    outfile<<"          "<<N_impropers_sp<<" impropers"<<endl;
    outfile<<endl;
    outfile<<"          "<<N_atom_types1<<" atom types"<<endl;
    outfile<<"          "<<N_bond_types1<<" bond types"<<endl;
    outfile<<"          "<<N_angle_types1<<" angle types"<<endl;
    outfile<<endl;
     
    if(fCombine){
      outfile<<"     "<<xb2lo<<"    "<<xb2hi<<"  xlo xhi"<<endl;
      outfile<<"     "<<yb2lo<<"    "<<yb2hi<<"  ylo yhi"<<endl;
      
      if(fShift)
	outfile<<"     "<<zb2lo<<"    "<<zb2hi + (zb1hi - zb1lo) <<"  zlo zhi"<<endl;
      else
	outfile<<"     "<<zb2lo<<"    "<<zb2hi<<"  zlo zhi"<<endl;

    } else {
//       outfile<<"     "<<xb1lo-xb1c<<"    "<<xb1hi-xb1c<<"  xlo xhi"<<endl;
//       outfile<<"     "<<yb1lo-yb1c<<"    "<<yb1hi-yb1c<<"  ylo yhi"<<endl;
      outfile<<"     "<<xb1lo<<"    "<<xb1hi<<"  xlo xhi"<<endl;
      outfile<<"     "<<yb1lo<<"    "<<yb1hi<<"  ylo yhi"<<endl;
      
      if(fShift)
	outfile<<"     "<<zb1lo<<"    "<<zb1hi<<"  zlo zhi"<<endl;
      else
	outfile<<"     "<<zb1lo<<"    "<<zb1hi<<"  zlo zhi"<<endl;
   
    }

    outfile<<endl;  
    outfile<<"Masses"<<endl;
    outfile<<endl;
    outfile<<"\t"<<1<<"\t"<<Mass[0]<<endl;
    outfile<<"\t"<<2<<"\t"<<Mass[1]<<endl;
    
 
  }
  
  //##############################
  //## generate atom list
  //##############################
  if(CheckDroplet){
    cout<<endl;
    cout<<"Atoms"<<endl;
    cout<<endl;
  }
  if(fWrite){
    outfile<<endl;
    outfile<<"Atoms"<<endl;
    outfile<<endl;
  }

  int nid,mid = 0,type;
  double ch,spx,spy,spz;

  for(int m = 0; m < N_atoms_sp; m++){
    
    //nid  = Atom_spid[m];
    //mid  = Atom_spmol[m];
    nid = m+1;
    type = Atom_sptype[m];
    spx  = Atom_spx[m];   
    spy  = Atom_spy[m];
    spz  = Atom_spz[m];
   
    if(m%N_atoms_mol1 == 0){
      mid++;

    }

    double dx = spx - xb1c;
    double dy = spy - yb1c;
    double dz;
    if(fFullSphere || fFullCylinder || fSlab)
      dz = spz - zb1c;
    else
      dz = spz;

    if(fShiftCentTozero){
      spx = dx;
      spy = dy;
      spz = dz;
    }
    //if(sqrt(dx*dx + dy*dy + dz*dz) > zb1c - zb1lo){
//     if(fabs(spz) > Rmax){
//       cout<<"spz is larger than Rmax for: "<<nid<<"\t"<<spx<<"\t"<<spy<<"\t"<<spz<<endl;
//     }
    if(m%3 == 0)
      ch = charge[0];
    else
      ch = charge[1];
    
    if(CheckDroplet){
      cout<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	  <<spx<<"  "<<spy<<"  "<<spz<<"  0  0  0"<<endl;	
    } 
    outfile<<"    "<<nid<<"    "<<mid<<"  "<<type<<"  "<<ch<<" "
	   <<spx<<"  "<<spy<<"  "<<spz<<endl;	
  
//       1    1  1 -0.8472   12.12456   28.09298   22.27452  0  1  0
//       2    1  2  0.4236   12.53683   28.75606   22.89928  0  1  0
//       3    1  2  0.4236   11.49482   28.56390   21.65678  0  1  0
  }

// //########################
// //## generate bond list
// //########################
  
  if(CheckDroplet){
    cout<<endl;
    cout<<"Bonds"<<endl;
    cout<<endl;
  }

  if(fWrite){
    outfile<<endl;
    outfile<<"Bonds"<<endl;
    outfile<<endl;
  }

   double cb = 1;

   for(int i = 0; i < N_bonds_sp; i++){

     if(CheckDroplet){
       if(i%2 == 0){
	 cout<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
       } else {
	 cout<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
       }
     }

    if(fWrite){
      if(i%2 == 0){
	outfile<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+1<<endl;
      } else {
	outfile<<"  "<<i+1<<"  "<<1<<"  "<<cb<<"  "<<cb+2<<endl;
      }
    }
    if(i%2 != 0)
      cb += 3;
    
  }

//       1   1      1      2
//       2   1      1      3
//       3   1      4      5
//       4   1      4      6

//#######################
//## generate angle list
//#######################

   if(CheckDroplet){
     cout<<endl;
     cout<<"Angles"<<endl;
     cout<<endl;
   }

  if(fWrite){
    outfile<<endl;
    outfile<<"Angles"<<endl;
    outfile<<endl;
  }

  double ca = 1;
 
  for(int i = 0; i < N_angles_sp; i++){

    if(CheckDroplet)
      cout<<"  "<<i+1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;

    if(fWrite)
      outfile<<"  "<<i+1<<"  "<<1<<"  "<<ca+1<<"  "<<ca<<"  "<<ca+2<<endl;
    
    ca += 3;
    
  }


//       1   1      2      1      3
//       2   1      5      4      6
//       3   1      8      7      9
//       4   1     11     10     12
//       5   1     14     13     15
//       6   1     17     16     18

   outfile.close();



}

void MakeWaterSphere(){

  fFullSphere = 1;
  fShiftCentTozero = 1;

  cout<<"########################################"<<endl;
  if(fFullSphere)
    cout<<"Cutting a sphere ..."<<endl;
  else
    cout<<"Cutting a hemisphere "<<endl;
  cout<<"########################################"<<endl;
  CalcBoxCenter();

  string name;


  if(fFullSphere)
    name = Form("waterdroplet.spce.%d",(Int_t)Rmax);
  else
    name = "waterhemisphere.spce.50";
    
  if(fShiftCentTozero){
    //Rmax = 50.;
    if(fFullSphere)
      name = Form("waterdroplet.shifted.to.zero.spce.%d",(Int_t)Rmax); 
    else
      name = "waterhemisphere.shifted.to.zero.spce.50";
    
    //ShiftCenterToZero();
  }
  if(fFullSphere)
    CutSphere(zb1c, n_time_steps1 - 1);
  else
    CutHemiSphere(zb1c, n_time_steps1 - 1);

  cout<<"########################################"<<endl;
  if(fFullSphere)
    cout<<"Done cutting a sphere "<<endl;
  else
    cout<<"Done cutting a hemisphere "<<endl;
  cout<<"########################################"<<endl;

  fWrite = 1;

  MakeSphereList(name);
  if(fCombine)
    CombineList(n_time_steps1 - 1);

}

void MakeWaterCylinder(){

  //fFullCylinder = 1;
  fShiftCentTozero = 1;

  cout<<"########################################"<<endl;
  if(fFullCylinder)
    cout<<"Cutting a cylinder ..."<<endl;
  else 
    cout<<"Cutting a Hemi-Cylinder "<<endl;
  cout<<"########################################"<<endl;
  CalcBoxCenter();

  TString name;

  if(fFullCylinder)
    name = Form("watercylinder.spce_yLo0.9_yHi49.1_%dA",(int)Rmax);
  else
    name = Form("waterhemicylinder.spce_yLo0.9_yHi49.1_%dA",(int)Rmax);
    
   if(fShiftCentTozero){
//     if(fFullCylinder)
//       name = Form("watercylinder.shifted.to.zero.spce_yLo0.9_yHi49.1_%dA",(int)Rmax); 
//     else
//       name = Form("waterhemicylinder.shifted.to.zero.spce_yLo0.9_yHi49.1_%dA",(int)Rmax);

    if(fFullCylinder)
       name = Form("watercylinder.shifted.to.zero.spce_%dA_y100",(int)Rmax);
    else
      name = Form("waterhemicylinder.shifted.to.zero.spce_%dA_y100",(int)Rmax);
    
     //ShiftCenterToZero();
  }
  if(fFullCylinder)
    CutCylinder(zb1c, n_time_steps1 - 1);
  else
    CutHemiCylinder(zb1c, n_time_steps1 - 1);

  cout<<"########################################"<<endl;
  if(fFullCylinder)
    cout<<"Done cutting a Cylinder "<<endl;
  else 
    cout<<"Done cutting a Hemi-Cylinder "<<endl;
  cout<<"########################################"<<endl;

  fWrite = 1;

  MakeSphereList(name.Data());
//   if(fCombine)
//     CombineList(n_time_steps1 - 1);

}

void MakeWaterSlab(){
  fSlab  = 1;
  fShiftCentTozero = 1;

  //xSlabLo += 1.0;
  xSlabHi -= 1.0;
  //ySlabLo += 1.0;
  ySlabHi -= 1.0;
  //zSlabLo;
  //zSlabHi;
  cout<<"########################################"<<endl;
  cout<<"Cutting a slab ..."<<endl;
  cout<<"########################################"<<endl;
  CalcBoxCenter();

  TString name;

  //name = Form("waterslab.shifted.to.zero.spce_x%dA_y%dA_z%dA",(int)(2*xSlabHi),(int)(2*ySlabHi),(int)(2*ySlabHi));
  //name = Form("waterslab3_x%dA_y%dA_z%dA.dat",(int)(2*xSlabHi),(int)(2*ySlabHi),(int)(2*ySlabHi));
  //name = Form("waterslab_sapphire_x50A_y50A_z100A.dat");
  //name = Form("waterslab_PMMA_x71A_y75A_z25A.dat");
  //name = Form("waterslab_PMMA_x2.5_73A_y-0.5_75.5A_z0_25A.dat");
  //name = Form("waterslab_PMMA_x0_75A_y0_75A_z0_25A.dat");
  name = Form("waterslab_for_CTAB_x0_35.36A_y0_35.36A_z0_27.52A.dat");

  //name = Form("CH4onTopOfH2Ofilm_Sel.dat");
      
  CutSlab(zb1c, n_time_steps1 - 1);

  cout<<"########################################"<<endl;
  cout<<"Done cutting a Slab "<<endl;
  cout<<"########################################"<<endl;

  fWrite = 1;

  MakeSphereList(name.Data());

}

void waterdroplet_tip4p_new(int rad = 40, TString obj = "fullCylinder"){

  if(obj.Contains("fullcylinder"))
     fFullCylinder = 1;
  if(obj.Contains("hemicylinder"))
     fFullCylinder = 0;
  if(obj.Contains("fullsphere"))
     fFullSphere = 1;
  if(obj.Contains("semisphere"))
     fFullSphere = 0;
  
  //usage:   ./a.out isdump radius
  //isdump: 0 for init file
  //isdump: 1 for dumpfle
  //radius: 42.5, 52.5 (in A)

  char buf[250];
  string index;

  Rmax = rad;   //in A
 
  //npt
//   int run = 1000000;
//   int dump_steps = 10000;
  //nvt

  //first file
  run1 = 1;//100000;
  dump_steps1 = 2;//10000;
  n_time_steps1 = 1;// + (run1/dump_steps1);

  //second file
  run2 = 1;//100000;
  dump_steps2 = 2;//10000;
  n_time_steps2 = 1;// + (run2/dump_steps2);
  
  Initialize();

  ReadFile1();
  //ReadFile2();

  DuplicateData();

  if(obj.Contains("sphere"))
     MakeWaterSphere();
  else if(obj.Contains("cylinder"))
    MakeWaterCylinder();
  else
    MakeWaterSlab();

}

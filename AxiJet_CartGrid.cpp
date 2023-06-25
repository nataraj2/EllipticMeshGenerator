#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
using namespace std;

const double pi = M_PI;

class Grid
{
	public:
		int Nx, Ny, Nz;
		double*** x;
		double*** y;
		double*** z;
		int*** iblank;
		void Iblanking(int nLev);

		double** rho	;
		double** ux;
		double** uy;
		double** p;
		double** E;
		double*** U;

		void set_grids(int n,int Npts_x,int Npts_y,const double xmin[],const double xmax[],const double ymin[],const double ymax[]);
};

double** Create2DMatrix(int, int);
int** Create2DMatrix_INT(int, int);
double*** Create3DMatrix(int, int, int);
int*** Create3DMatrix_INT(int, int, int);
double**** Create4DMatrix(int, int, int, int);
void Delete4DMatrix(double****&, int, int, int, int);
void Delete3DMatrix(double***&, int, int, int);
void Delete3DMatrix_INT(int***&, int, int, int);
void Delete2DMatrix(double**&, int, int);

double compute_coord_val(int i, int nx, double delta_x, double pos_x, double coord_min, double coord_max);
void create_AxiJet_CartGrid(int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax);
void write_plot3d_grid(double*** x,double*** y,double*** z,int Nx,int Ny,int Nz, std::string filename_string);


int main()
{
	int Nblocks = 3;
	Grid* grid = new Grid [Nblocks];

	double xmin = 0.0;
	double xmax = 15.0;
	
	double ymin = 0.0;
	double ymax = 15.0;

	int nx = 300, ny = 200, nz = 1;
	
	create_AxiJet_CartGrid(nx, ny, nz, xmin, xmax, ymin, ymax);

	return 0;
}

double compute_coord_val(int i, int nx, double delta_x, double pos_x, double coord_min, double coord_max)
{
	pos_x = pos_x/(coord_max-coord_min);
	
	double xi = float(i)/float(nx-1);
	double u1 = tanh(delta_x*(1.0-xi))*(1.0-pos_x);
	double u2 = (2.0-tanh(delta_x*xi))*pos_x;
	double fac = 1.0 - ((u1+u2)-pos_x);

	double coord_val = coord_min + fac*(coord_max-coord_min);

	return coord_val;

}

void create_AxiJet_CartGrid(int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax)
{

	int npoints_imp_wall_zone = 100;

	double*** x = Create3DMatrix(nx+npoints_imp_wall_zone,ny,nz);
    double*** y = Create3DMatrix(nx+npoints_imp_wall_zone,ny,nz);
    double*** z = Create3DMatrix(nx+npoints_imp_wall_zone,ny,nz);

	int ipos2;
	double xpos2 = 11.0;
	bool found_index = false;


	double delta_y = 6.0;
	double pos_y = 1.0;

	for(int i=0;i<nx; i++){

		double xval = compute_coord_val(i,nx,3.0,6.0,xmin,xmax);
	
		if(xval > xpos2 and !found_index){
			ipos2 = i;
			found_index = true;
			break;
		}
			
		for(int j=0;j<ny;j++){
			double yval = compute_coord_val(j,ny,delta_y,pos_y,ymin,ymax); 
			for(int k=0;k<nz;k++){
				x[i][j][k] = xval;
				y[i][j][k] = yval;
				z[i][j][k] = 0.0;
			}
		}
	}
	
	for(int i=ipos2;i<nx+npoints_imp_wall_zone; i++){
	
		double xval = compute_coord_val(i-ipos2,nx+npoints_imp_wall_zone-ipos2,2.0,3.0,xpos2,xmax);

		for(int j=0;j<ny;j++){
			double yval = compute_coord_val(j,ny,delta_y,pos_y,ymin,ymax); 
			for(int k=0;k<nz;k++){
				x[i][j][k] = xval;
				y[i][j][k] = yval;
				z[i][j][k] = 0.0;
			}
		}
	}
	
	write_plot3d_grid(x,y,z,nx+npoints_imp_wall_zone,ny,nz,"axijet_mesh.xyz");
		
}

void write_plot3d_grid(double*** x,double*** y,double*** z,int Nx,int Ny,int Nz, std::string filename_string)
{
	ofstream outFile;
	const char* filename = filename_string.c_str();
	outFile.open (filename, ios::out | ios::binary);

	int size = 4;
	int Nblocks = 1;
	outFile.write((char*) &size, sizeof(int));
	outFile.write((char*) &Nblocks, sizeof(Nblocks));
	outFile.write((char*) &size, sizeof(int));

	size = 4*Nblocks*3;
	outFile.write((char*) &size, sizeof(int));
	for (int i=0;i<Nblocks;i++)
	{			
		outFile.write((char*) &Nx, sizeof(Nx));
		outFile.write((char*) &Ny, sizeof(Ny));
		outFile.write((char*) &Nz, sizeof(Nz));
	}
	outFile.write((char*) &size, sizeof(int));

	for (int ii=0;ii<Nblocks;++ii)
	{
		size = 8*3*Nx*Ny*Nz + 4*Nx*Ny*Nz;
		outFile.write((char*) &size, sizeof(int));
		for (int k=0;k<Nz;++k)
			for (int j=0;j<Ny;++j)
				for (int i=0;i<Nx;++i)
				{
					outFile.write((char*) &x[i][j][k], sizeof(x[i][j][k]));
				}

		for (int k=0;k<Nz;++k)
			for (int j=0;j<Ny;++j)
				for (int i=0;i<Nx;++i)
				{
					outFile.write((char*) &y[i][j][k], sizeof(y[i][j][k]));
				}

		for (int k=0;k<Nz;++k)
			for (int j=0;j<Ny;++j)
				for (int i=0;i<Nx;++i)
				{
					outFile.write((char*) &z[i][j][k], sizeof(z[i][j][k]));
				}

		for (int k=0;k<Nz;++k)
			for (int j=0;j<Ny;++j)
				for (int i=0;i<Nx;++i)
				{
					int tmp = 1;
					outFile.write((char*) &tmp, sizeof(int));
				}

		outFile.write((char*) &size, sizeof(int));
	}

	outFile.close();
	cout << "===============================" << endl;
	printf("Wrote \"%s\" \n", filename);
	cout << "===============================" << endl;
}

double** Create2DMatrix(int Ni, int Nj)
{
    double **the_array = new double* [Ni];
    double *tempxy = new double[Ni*Nj];
 		for ( int i = 0 ; i < Ni; ++i, tempxy += Nj ) {
			the_array[i] = tempxy;
    }

    for(int i(0); i < Ni; ++i)
        for(int j(0); j < Nj; ++j)
		{	the_array[i][j]= 0.;		}

    /*double** the_array = new double* [Ni];
    for(int i(0); i < Ni; ++i)
    {
        the_array[i] = new double[Nj];

        for(int j(0); j < Nj; ++j)
        {
            the_array[i][j] = 0;            
        }
    }*/

    return the_array;
}

int** Create2DMatrix_INT(int Ni, int Nj)
{
    int **the_array = new int* [Ni];
    int *tempxy = new int[Ni*Nj];
 		for ( int i = 0 ; i < Ni; ++i, tempxy += Nj ) {
			the_array[i] = tempxy;
    }

    for(int i(0); i < Ni; ++i)
        for(int j(0); j < Nj; ++j)
		{	the_array[i][j]= 0;		}

    /*int** the_array = new int* [Ni];
    for(int i(0); i < Ni; ++i)
    {
        the_array[i] = new int[Nj];

        for(int j(0); j < Nj; ++j)
        {
            the_array[i][j] = 0;            
        }
    }*/

    return the_array;
}

double*** Create3DMatrix(int Nk, int Ni, int Nj)
{
    double ***the_array = new double**[Nk];
    double **tempxy = new double*[Nk*Ni];
    double *tempxyz = new double[Nk*Ni*Nj];
    for ( int k = 0 ; k < Nk ; ++k, tempxy += Ni ) {
        the_array[k] = tempxy;
		for ( int i = 0 ; i < Ni; ++i, tempxyz += Nj ) {
			the_array[k][i] = tempxyz;
    } }

    for(int k(0); k < Nk; ++k)
        for(int i(0); i < Ni; ++i)
            for(int j(0); j < Nj; ++j)
			{	the_array[k][i][j]= 0.;		}

    return the_array;
}

double**** Create4DMatrix(int Nl, int Ni, int Nj, int Nk)
{
    double ****the_array = new double***[Nl];
    double ***tempxy = new double**[Nl*Ni];
    double **tempxyz = new double*[Nl*Ni*Nj];
    double *tempxyzl = new double[Nl*Ni*Nj*Nk];
    for ( int l = 0 ; l < Nl ; ++l, tempxy += Ni ) {
        the_array[l] = tempxy;
		for ( int i = 0 ; i < Ni; ++i, tempxyz += Nj ) {
			the_array[l][i] = tempxyz;
			for ( int j = 0 ; j < Nj; ++j, tempxyzl += Nk ) {
				the_array[l][i][j] = tempxyzl;
    } } }

    for(int l(0); l < Nl; ++l)
	    for(int i(0); i < Ni; ++i)
	        for(int j(0); j < Nj; ++j)
			    for(int k(0); k < Nk; ++k)
				{	the_array[l][i][j][k]= 0.;		}

    return the_array;
}

int*** Create3DMatrix_INT(int Nk, int Ni, int Nj)
{
    int ***the_array = new int**[Nk];
    int **tempxy = new int*[Nk*Ni];
    int *tempxyz = new int[Nk*Ni*Nj];
    for ( int k = 0 ; k < Nk ; ++k, tempxy += Ni ) {
        the_array[k] = tempxy;
		for ( int i = 0 ; i < Ni; ++i, tempxyz += Nj ) {
			the_array[k][i] = tempxyz;
    } }

    for(int k(0); k < Nk; ++k)
        for(int i(0); i < Ni; ++i)
            for(int j(0); j < Nj; ++j)
			{	the_array[k][i][j]= 0.;		}

    return the_array;
}

void Delete4DMatrix(double****& the_array, int Nl, int Ni, int Nj, int Nk)
{
    delete [] the_array[0][0][0];
    delete [] the_array[0][0];
    delete [] the_array[0];
    delete [] the_array;
}

void Delete3DMatrix(double***& the_array, int Nk, int Ni, int Nj)
{
    delete [] the_array[0][0];
    delete [] the_array[0];
    delete [] the_array;
}

void Delete3DMatrix_INT(int***& the_array, int Nk, int Ni, int Nj)
{
    delete [] the_array[0][0];
    delete [] the_array[0];
    delete [] the_array;
}

void Delete2DMatrix(double**& the_array, int Ni, int Nj)
{
    delete [] the_array[0];
    delete [] the_array;
}

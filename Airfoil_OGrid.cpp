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

void get_airfoil_coords(double*&,double*&,int);
void write_plot3d_grid(double*** x,double*** y,double*** z,int Nx,int Ny,int Nz, std::string filename_string);
void write_plot3d_BLOCK_GRID(Grid* grid, int Nblocks);
void create_mesh(int nx, int ny, int nz, double***& x, double***&y, double***& z,double zmin, double zmax, 
				 const double *x_xi_0, const double *x_xi_max, const double *x_eta_0, const double *x_eta_max,
				 const double *y_xi_0, const double *y_xi_max, const double *y_eta_0, const double *y_eta_max, double delta_y, 
				 int n_iterations, bool invert, std::string filename);

void create_OGrid_mesh(Grid* gl, int nx, int ny, int nz, double zmin, double zmax, int indx[2], double* x_afoil, double* y_afoil, 
					   double flat_portion, double semicircle_portion, double x_circle_start, double len_y, double delta_y, int n_iterations);

void create_top_right_mesh(Grid* gl, int nx, int ny, int nz, int n_add, double zmin, double zmax, int indx[2], 
							double* x_afoil, double* y_afoil, double len_x, double len_y, double delta_y, int n_iterations);


void create_bottom_right_mesh(Grid* gl, int nx, int ny, int nz, int n_add, double zmin, double zmax, int indx[2], 
							  double* x_afoil, double* y_afoil, double len_x, double len_y, int npts_afoil, double delta_y, int n_iterations);

int main()
{
	int Nblocks = 3;
	Grid* grid = new Grid [Nblocks];
	
	double zmin = 0.0, zmax = 1.0;
	int nz = 20;

	int npts_afoil = 100;			
	int ny = 100;

	// Number of smoothing iterations for the elliptic smoothing
	int n_iterations = 20;

	// Choose the percentage of chord from the leading edge for the 
	// C-grid portion

	double c_grid_portion = 100.0;

	// Parameters for the O grid

	double flat_portion = 25.0; // Should be less than 50%  - the total percent of points on the flat portion
	double semicircle_portion = 25.0;
	
	if(2*flat_portion + semicircle_portion > 100.0){
		cout << "Twice flat portion plus semicircle portion is " << 2*flat_portion + semicircle_portion << 
				" but cannot be greater than or equal to 100.0...Exiting" << "\n";
		exit(0);
	}
		
	double len_y = 1.0;		// radius of the circle
	double x_circle_start = 0.0;	// center of the circle

	// Additional length in streamwise direction beyond the trailing edge
	double len_x = 1.0;		// only used on the top and the bottom right grid

	double delta_y = 2.0;	// stretching in the y direction

	// Parameter list ends

	double* x_afoil	= new double[npts_afoil];
	double* y_afoil	= new double[npts_afoil];
	get_airfoil_coords(x_afoil,y_afoil,npts_afoil);

		
	// Find the co-ordinates of the airfoil on the top and bottom which correspond to this value

	double val = c_grid_portion/100.0*1.0;
	int indx[2];
	int count = 0;
	for(int i=0;i<npts_afoil-1;i++){
		if((x_afoil[i]-val)*(x_afoil[i+1]-val) < 0){
			//indx[count] = i;
			//count++;
		}
	}

	indx[0] = 0;
	indx[1] = npts_afoil-1;
	int nx = indx[1] - indx[0] + 1;		// for the C grid

	// Create the C grid

	create_OGrid_mesh(&grid[0], nx, ny, nz, zmin, zmax, indx, x_afoil, y_afoil, flat_portion, semicircle_portion, x_circle_start, len_y, delta_y, n_iterations);

	// Create the top right mesh
	
	// Find the coordinates of the airfoil on the top and bottom which correspond to this value
	double pct_afoil = 70.0;
	val = pct_afoil/100.0*1.0;
	count = 0;
	for(int i=0;i<npts_afoil-1;i++){
		if((x_afoil[i]-val)*(x_afoil[i+1]-val) < 0){
			indx[count] = i;
			count++;
		}
	}

	indx[0] = indx[0]+1;

	std::cout << "The airfoil surface is from 1 to " << indx[0]+1 << "\n";

	int n_points_flat = flat_portion/100.0*nx;
	double delx_flat = (x_afoil[indx[0]]-x_circle_start)/(n_points_flat-1);
	
	int total_nx = (1.0 + (1.0 - x_afoil[indx[0]]))/delx_flat;

	int n_add = total_nx - indx[0];
	nx = total_nx; // Additional points

	create_top_right_mesh(&grid[1], nx, ny, nz, n_add, zmin, zmax, indx, x_afoil, y_afoil, len_x, len_y, delta_y, n_iterations);

	// Create the bottom right mesh

	create_bottom_right_mesh(&grid[2], nx, ny, nz, n_add, zmin, zmax, indx, x_afoil, y_afoil, len_x, len_y, npts_afoil, delta_y, n_iterations);
	
	write_plot3d_BLOCK_GRID(grid,Nblocks);

	return 0;
}

void write_plot3d_BLOCK_GRID(Grid* grid, int Nblocks)
{
	ofstream outFile;
	char filename[] = "afoil_block_grid.xyz";
	outFile.open (filename, ios::out | ios::binary);

	int size = 4;
	outFile.write((char*) &size, sizeof(int));
	outFile.write((char*) &Nblocks, sizeof(Nblocks));
	outFile.write((char*) &size, sizeof(int));

	size = 4*Nblocks*3;
	outFile.write((char*) &size, sizeof(int));
	for (int i=0;i<Nblocks;i++)
	{			
		outFile.write((char*) &grid[i].Nx, sizeof(grid[i].Nx));
		outFile.write((char*) &grid[i].Ny, sizeof(grid[i].Ny));
		outFile.write((char*) &grid[i].Nz, sizeof(grid[i].Nz));
	}
	outFile.write((char*) &size, sizeof(int));

	for (int ii=0;ii<Nblocks;++ii)
	{
		size = 8*3*grid[ii].Nx*grid[ii].Ny*grid[ii].Nz + 4*grid[ii].Nx*grid[ii].Ny*grid[ii].Nz;
		outFile.write((char*) &size, sizeof(int));
		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].x[i][j][k], sizeof(grid[ii].x[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].y[i][j][k], sizeof(grid[ii].y[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].z[i][j][k], sizeof(grid[ii].z[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
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

void create_bottom_right_mesh(Grid* gl, int nx, int ny, int nz, int n_add, double zmin, double zmax, int indx[2], 
							double* x_afoil, double* y_afoil, double len_x, double len_y, int npts_afoil, double delta_y, int n_iterations)
{
	double*** x = Create3DMatrix(nx,ny,nz);
	double*** y = Create3DMatrix(nx,ny,nz);
	double*** z = Create3DMatrix(nx,ny,nz);

	// Create the boundary values for x and y

	double* x_xi_0 	   = new double[ny];
	double* x_xi_max   = new double[ny];
	double* x_eta_0    = new double[nx];
	double* x_eta_max  = new double[nx];
	double* y_xi_0     = new double[ny];
	double* y_xi_max   = new double[ny]; 
	double* y_eta_0    = new double[nx]; 
	double* y_eta_max  = new double[nx];
	
	// Define eta_0 values
	// There are indx[0] + 1 points on the airfoil surface on the top


	double delx_bottom = len_x/(n_add-1); 
	double delx_top = (len_x + (1.0 - x_afoil[indx[1]]))/(nx-1);

	for(int i = 0; i<nx; i++){
		if(indx[1]+i<=npts_afoil-1){
			x_eta_max[i] = x_afoil[indx[1]+i];
			y_eta_max[i] = y_afoil[indx[1]+i];
		}
		else{
			x_eta_max[i] = 1.0 + (i-indx[0])*delx_bottom;
			y_eta_max[i] = 0.0;
		}
		x_eta_0[i] = x_afoil[indx[0]] + i*delx_top;
		y_eta_0[i] = -len_y; 	
	}


	
	for(int j=0;j<ny; j++){
		double eta = float(j)/float(ny-1);
		double fac = 1.0 + tanh(delta_y*(eta-1))/tanh(delta_y);
		y_xi_0[ny-1-j] = y_afoil[indx[1]] - fac*(len_y + y_afoil[indx[1]]);
		y_xi_max[ny-1-j] = 0 + fac*(-len_y);

		x_xi_0[j] = x_afoil[indx[0]];
		x_xi_max[j] = 1.0 + len_x;
	}
	create_mesh(nx,ny,nz,x,y,z,zmin,zmax,
				x_xi_0,x_xi_max,x_eta_0,x_eta_max,
				y_xi_0,y_xi_max,y_eta_0,y_eta_max, delta_y,  
				n_iterations,true,"file_mesh_bottom_right.vtk");


	write_plot3d_grid(x,y,z,nx,ny,nz,"mesh_bottom.xyz");
	
	gl->Nx = nx;
	gl->Ny = ny;
	gl->Nz = nz;
	gl->x = Create3DMatrix(nx,ny,nz);
	gl->y = Create3DMatrix(nx,ny,nz);
	gl->z = Create3DMatrix(nx,ny,nz);
	for (int k=0;k<nz;++k)
		for (int j=0;j<ny;++j)
			for (int i=0;i<nx;++i)
			{
				gl->x[i][j][k] = x[i][j][k];
				gl->y[i][j][k] = y[i][j][k];
				gl->z[i][j][k] = z[i][j][k];
			}
			
	Delete3DMatrix(x, nx, ny, nz);
	Delete3DMatrix(y, nx, ny, nz);
	Delete3DMatrix(z, nx, ny, nz);
}

void create_top_right_mesh(Grid* gl, int nx, int ny, int nz, int n_add, double zmin, double zmax, int indx[2], 
						   double* x_afoil, double* y_afoil, double len_x, double len_y, double delta_y, int n_iterations)
{

	double*** x = Create3DMatrix(nx,ny,nz);
	double*** y = Create3DMatrix(nx,ny,nz);
	double*** z = Create3DMatrix(nx,ny,nz);

	// Create the boundary values for x and y

	double* x_xi_0 	   = new double[ny];
	double* x_xi_max   = new double[ny];
	double* x_eta_0    = new double[nx];
	double* x_eta_max  = new double[nx];
	double* y_xi_0     = new double[ny];
	double* y_xi_max   = new double[ny]; 
	double* y_eta_0    = new double[nx]; 
	double* y_eta_max  = new double[nx];
	
	// Define eta_0 values
	// There are indx[0] + 1 points on the airfoil surface on the top


	double delx_bottom = len_x/(n_add-1); 
	double delx_top = (len_x + (1.0 - x_afoil[indx[0]]))/(nx-1);

	for(int i = 0; i<nx; i++){
		if(indx[0]-i>=0){
			x_eta_0[i] = x_afoil[indx[0]-i];
			y_eta_0[i] = y_afoil[indx[0]-i];
		}
		else{
			x_eta_0[i] = 1.0 + (i-indx[0])*delx_bottom;
			y_eta_0[i] = 0.0;
		}
		x_eta_max[i] = x_afoil[indx[0]] + i*delx_top;
		y_eta_max[i] = len_y; 	
	}


	for(int j=0;j<ny; j++){
		double eta = float(j)/float(ny-1);
		double fac = 1.0 + tanh(delta_y*(eta-1))/tanh(delta_y);
		y_xi_0[j] = y_afoil[indx[0]] + fac*(len_y-y_afoil[indx[0]]);
		y_xi_max[j] = fac*(len_y);

		x_xi_0[j] = x_afoil[indx[0]];
		x_xi_max[j] = 1.0 + len_x;
	}
	create_mesh(nx,ny,nz,x,y,z,zmin,zmax,
				x_xi_0,x_xi_max,x_eta_0,x_eta_max,
				y_xi_0,y_xi_max,y_eta_0,y_eta_max, delta_y,  
				n_iterations,false,"file_mesh_top_right.vtk");


	write_plot3d_grid(x,y,z,nx,ny,nz,"mesh_top.xyz");
	
	gl->Nx = nx;
	gl->Ny = ny;
	gl->Nz = nz;
	gl->x = Create3DMatrix(nx,ny,nz);
	gl->y = Create3DMatrix(nx,ny,nz);
	gl->z = Create3DMatrix(nx,ny,nz);
	for (int k=0;k<nz;++k)
		for (int j=0;j<ny;++j)
			for (int i=0;i<nx;++i)
			{
				gl->x[i][j][k] = x[i][j][k];
				gl->y[i][j][k] = y[i][j][k];
				gl->z[i][j][k] = z[i][j][k];
			}
			
	Delete3DMatrix(x, nx, ny, nz);
	Delete3DMatrix(y, nx, ny, nz);
	Delete3DMatrix(z, nx, ny, nz);
}

void create_OGrid_mesh(Grid* gl, int nx, int ny, int nz, double zmin, double zmax, int indx[2], double* x_afoil, double* y_afoil, 
					   double flat_portion, double semicircle_portion, double x_circle_start, double len_y, double delta_y, int n_iterations)
{
	// Define the coordinate arrays
	double*** x = Create3DMatrix(nx,ny,nz);
	double*** y = Create3DMatrix(nx,ny,nz);
	double*** z = Create3DMatrix(nx,ny,nz);

	// Create the boundary values for x and y

	double* x_xi_0 	   = new double[ny];
	double* x_xi_max   = new double[ny];
	double* x_eta_0    = new double[nx];
	double* x_eta_max  = new double[nx];
	double* y_xi_0     = new double[ny];
	double* y_xi_max   = new double[ny]; 
	double* y_eta_0    = new double[nx]; 
	double* y_eta_max  = new double[nx];
	
	// Assign values to the boundary values of x and y

	// Define eta_0 values 
	for(int i=0;i<nx; i++){
		x_eta_0[i] = x_afoil[nx-1-i];
		y_eta_0[i] = y_afoil[nx-1-i];
 	}

	// Define eta_max values i.e the outer domain of the C grid 
	// Assign flat_line_portion, circle portion (there are 2 flat line portions)
;
	int n_points_flat = flat_portion/100.0*nx;

	// The flat portion extends from x_circle_start to x_eta_0[0] to


	// Circle

	double quartercircle_portion = (100.0 - 2*flat_portion - semicircle_portion)/2.0;
	int n_quartercircle = quartercircle_portion/100.0*nx;

	int n_semicircle = nx - 2*(n_points_flat+n_quartercircle-1)+2;

	for(int i=0;i<nx; i++){
	
		// Bottom quarter circle	
		if(i< n_quartercircle){
			double del_theta = (pi/2.0)/(n_quartercircle-1); 
			x_eta_max[i] = 1.0 + len_y*cos(2.0*pi - del_theta*i);
			y_eta_max[i] = len_y*sin(2.0*pi - del_theta*i);
		}
	
		// Bottom flat
		else if(i< n_quartercircle + n_points_flat - 1){
			int istart = n_quartercircle-1;
			double delx_flat = (1.0-x_circle_start)/(n_points_flat-1);		
			x_eta_max[i] = 1.0 - (i-istart)*delx_flat;
			y_eta_max[i] = -len_y;
		}
		// Semi circle
		else if (i< n_quartercircle + n_points_flat + n_semicircle - 2){
			int istart = n_quartercircle + n_points_flat-2;
			double del_theta = pi/(n_semicircle-1);
			x_eta_max[i] = x_circle_start + len_y*cos(3.0*pi/2.0 - del_theta*(i-istart));
			y_eta_max[i] = len_y*sin(3.0*pi/2.0 - del_theta*(i-istart));
		}
		// Top flat
		else if(i < n_quartercircle + n_points_flat + n_semicircle + n_points_flat -3){
			int istart = n_quartercircle + n_points_flat + n_semicircle - 3;
			double delx_flat = (1.0-x_circle_start)/(n_points_flat-1);
			x_eta_max[i] = x_circle_start + (i-istart)*delx_flat;
			y_eta_max[i] = len_y;
		}
		// Top quarter circle
		else{
			int istart = n_quartercircle + n_points_flat + n_semicircle + n_points_flat -4;
			double del_theta = (pi/2.0)/(n_quartercircle-1); 
			x_eta_max[i] = 1.0 + len_y*cos(pi/2 - del_theta*(i-istart));
			y_eta_max[i] = len_y*sin(pi/2 - del_theta*(i-istart));
		}
		

	}			

	// Assign xi_0, xi_max values of x and y

	//double dely_xi_0 = (y_afoil[indx[1]] - (-1.0*len_y))/(ny-1);
	//double dely_xi_max = (len_y - y_afoil[indx[0]])/(ny-1);

	for(int j=0;j<ny; j++){
		double eta = float(j)/float(ny-1);
		double fac = 1.0 + tanh(delta_y*(eta-1))/tanh(delta_y);
		y_xi_0[j] = 0.0;//y_afoil[indx[1]] + fac*(-len_y-y_afoil[indx[1]]);
		y_xi_max[j] = 0.0;//y_afoil[indx[0]] + fac*(len_y-y_afoil[indx[0]]);

		x_xi_0[j] = 1.0 + len_y/(ny-1)*j;
		x_xi_max[j] = 1.0 + len_y/(ny-1)*j;
	}

	// Create the mesh by passing the boundary definition (i.e. xi_0, xi_max, eta_0, eta_max) 

	create_mesh(nx,ny,nz,x,y,z,zmin,zmax,
				x_xi_0,x_xi_max,x_eta_0,x_eta_max,
				y_xi_0,y_xi_max,y_eta_0,y_eta_max, delta_y,  
				n_iterations,false,"file_mesh.vtk");

	write_plot3d_grid(x,y,z,nx,ny,nz,"mesh_OGrid.xyz");
	
	gl->Nx = nx;
	gl->Ny = ny;
	gl->Nz = nz;
	gl->x = Create3DMatrix(nx,ny,nz);
	gl->y = Create3DMatrix(nx,ny,nz);
	gl->z = Create3DMatrix(nx,ny,nz);
	for (int k=0;k<nz;++k)
		for (int j=0;j<ny;++j)
			for (int i=0;i<nx;++i)
			{
				gl->x[i][j][k] = x[i][j][k];
				gl->y[i][j][k] = y[i][j][k];
				gl->z[i][j][k] = z[i][j][k];
			}
			
	Delete3DMatrix(x, nx, ny, nz);
	Delete3DMatrix(y, nx, ny, nz);
	Delete3DMatrix(z, nx, ny, nz);
}

void create_mesh(int nx, int ny, int nz, double***& x, double***&y, double***& z, double zmin, double zmax, 
				 const double *x_xi_0, const double *x_xi_max, const double *x_eta_0, const double *x_eta_max,
				 const double *y_xi_0, const double *y_xi_max, const double *y_eta_0, const double *y_eta_max, double delta_y, 
				 int n_iterations, bool invert, std::string filename)
{

	FILE* file_mesh_vtk;
	file_mesh_vtk = fopen(filename.c_str(), "w");
	fprintf(file_mesh_vtk,"%s\n", "# vtk DataFile Version 2.0");
	fprintf(file_mesh_vtk,"%s\n", "Simple example");
	fprintf(file_mesh_vtk,"%s\n", "ASCII");
	fprintf(file_mesh_vtk,"%s\n", "DATASET STRUCTURED_GRID");
	fprintf(file_mesh_vtk,"%s %d %d %d\n", "DIMENSIONS", nx, ny, nz);
	fprintf(file_mesh_vtk,"%s %d %s\n", "POINTS", nx*ny*nz, "float");

	cout << "The mesh dimensions is x, y, z are " << nx << " " << ny << " " << nz << "\n";

	cout << "Initializing mesh generation" << "\n";	

	for(int k=0;k<nz;k++){
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
				x[i][j][k] = 0.0;
				y[i][j][k] = 0.0;
				z[i][j][k] = 0.0;
			}
		}
	}

	for(int k=0;k<nz;k++){
		for(int i=0;i<nx; i++){
	
			x[i][0][k] = x_eta_0[i];
			x[i][ny-1][k] = x_eta_max[i];
			y[i][0][k] = y_eta_0[i];
			y[i][ny-1][k] = y_eta_max[i];
	 	}

		for(int j=0;j<ny; j++){
			y[0][j][k] = y_xi_0[j];
			y[nx-1][j][k] = y_xi_max[j];
			x[0][j][k] = x_xi_0[j];
			x[nx-1][j][k] = x_xi_max[j];
		}
	}

	// Stretch in z
	double delta = 4.0;
	double actual_pos = 0.5;
	double pos = 1.0-actual_pos;

	for(int k=0;k<nz;k++){
		double zeta = float(k)/float(nz-1);
		double u1 = tanh(delta*(1.0-zeta))*(1.0-pos);
		double u2 = (2.0-tanh(delta*zeta))*pos;
		double fac = 1.0 - ((u1+u2)-pos);
		
		double delz = (zmax-zmin)/(nz-1);

		for(int i=0;i<nx; i++){
			for(int j=0;j<ny;j++){
				z[i][j][k] = zmin + k*delz;
			}
		}
	}	
	

	// Create the initial mesh slice-by-slice in the z direction
	for(int k=0;k<nz;k++){

		// Do transfinite interpolation to fill the mesh in the domain interior using the boundary values of x and y defined above
		int m = nx-1;
		int n = ny-1;

		for(int i=1;i<nx-1; i++){
			for(int j=1;j<ny-1;j++){
				x[i][j][k] = float(i)/float(m)*x[m][j][k] + float(m-i)/float(m)*x[0][j][k] + float(j)/float(n)*x[i][n][k] + float(n-j)/float(n)*x[i][0][k] - float(i)/float(m)*float(j)/float(n)*x[m][n][k] - 
					  float(i)/float(m)*float(n-j)/float(n)*x[m][0][k] - float(m-i)/float(m)*float(j)/float(n)*x[0][n][k] - float(m-i)/float(m)*float(n-j)/float(n)*x[0][0][k];

			//	y[i][j][k] = float(i)/float(m)*y[m][j][k] + float(m-i)/float(m)*y[0][j][k] + float(j)/float(n)*y[i][n][k] + float(n-j)/float(n)*y[i][0][k] - float(i)/float(m)*float(j)/float(n)*y[m][n][k] - 
			//		  float(i)/float(m)*float(n-j)/float(n)*y[m][0][k] - float(m-i)/float(m)*float(j)/float(n)*y[0][n][k] - float(m-i)/float(m)*float(n-j)/float(n)*y[0][0][k];

				double eta = float(j)/float(ny-1);
				double fac = 1.0 + tanh(delta_y*(eta-1))/tanh(delta_y);
				if(!invert){
					y[i][j][k] = y[i][0][k] + fac*(y[i][n][k]-y[i][0][k]);
				}
				else{
					y[i][n-j][k] = y[i][n][k] + fac*(y[i][0][k]-y[i][n][k]);;
				}
				
			}
		}
	
	} // End of k loop. Initial mesh generation for the whole domain ends here for the entire 3d domain

	cout << "Done with initial mesh generation" << "\n";


	// Do elliptic smoothing. Smooth only y. Elliptic equation is solved using Gauss-Siedel iteration
	for(int iter=0;iter<n_iterations;iter++){ 

		cout << "Doing smoothing iteration " << iter+1 << "\n";
		for(int i=1;i<nx-1;i++){
			for(int j=1;j<ny-1;j++){
				for(int k=0;k<nz;k++){

					double dxdxi = (x[i+1][j][k]-x[i-1][j][k])/2.0;
					double dydxi = (y[i+1][j][k]-y[i-1][j][k])/2.0;
					double dzdxi = (z[i+1][j][k]-z[i-1][j][k])/2.0;
		

					double dxdeta = (x[i][j+1][k]-x[i][j-1][k])/2.0;
					double dydeta = (y[i][j+1][k]-y[i][j-1][k])/2.0;
					double dzdeta = (z[i][j+1][k]-z[i][j-1][k])/2.0;

					double dxdzeta = (x[i][j][k+1]-x[i][j][k-1])/2.0;
					double dydzeta = (y[i][j][k+1]-y[i][j][k-1])/2.0;
					double dzdzeta = (z[i][j][k+1]-z[i][j][k-1])/2.0;

					double a_11, a_12, a_13, a_22, a_23, a_33;
					double alpha_11, alpha_12, alpha_13, alpha_22, alpha_23, alpha_33;

					a_11 = dxdxi*dxdxi    + dydxi*dydxi ;
					a_12 = dxdxi*dxdeta   + dydxi*dydeta ;
					a_13 = 0.0;
					a_22 = dxdeta*dxdeta   + dydeta*dydeta;
					a_23 = 0.0 ;
					a_33 = 0.0;

					alpha_11 = a_22;
					alpha_12 = -a_12;
					alpha_13 = 0.0;
					alpha_22 = a_11;
					alpha_23 = 0.0;
					alpha_33 = 0.0;
	
					double d2xdxi2_part  = (x[i+1][j][k]+x[i-1][j][k]);
					double d2ydxi2_part  = (y[i+1][j][k]+y[i-1][j][k]);
					double d2zdxi2_part  = (z[i+1][j][k]+z[i-1][j][k]);
				
					double d2xdeta2_part = (x[i][j+1][k]+x[i][j-1][k]);
					double d2ydeta2_part = (y[i][j+1][k]+y[i][j-1][k]);
					double d2zdeta2_part = (z[i][j+1][k]+z[i][j-1][k]);

					double d2xdzeta2_part  = (x[i][j][k+1]+x[i][j][k-1]);
					double d2ydzeta2_part  = (y[i][j][k+1]+y[i][j][k-1]);
					double d2zdzeta2_part  = (z[i][j][k+1]+z[i][j][k-1]);
	

					double d2xdxideta  = ((x[i+1][j+1][k]-x[i-1][j+1][k])/2.0 - (x[i+1][j-1][k]-x[i-1][j-1][k])/2.0)/2.0;

					double d2ydxideta =  ((y[i+1][j+1][k]-y[i-1][j+1][k])/2.0 - (y[i+1][j-1][k]-y[i-1][j-1][k])/2.0)/2.0;
					double d2ydxidzeta = ((y[i+1][j][k+1]-y[i-1][j][k+1])/2.0 - (y[i+1][j][k-1]-y[i-1][j][k-1])/2.0)/2.0; 
					
	
					double fac = 2.0*(alpha_11+alpha_22);
					
		
					y[i][j][k] = (alpha_11*d2ydxi2_part + 2.0*alpha_12*d2ydxideta + alpha_22*d2ydeta2_part + alpha_33*d2ydzeta2_part )/fac;
					x[i][j][k] = (alpha_11*d2xdxi2_part + 2.0*alpha_12*d2xdxideta + alpha_22*d2xdeta2_part)/fac;

				}
			}
		}
	}
			
	// Write VTK file
	//for(int k=0; k<nz; k++){
    //   	for(int j=0; j<ny; j++){
	//		for(int i=0; i<nx; i++){
	//			fprintf(file_mesh_vtk,"%g %g %g\n", x[i][j][k], y[i][j][k], z[i][j][k]);
	//		}
	//	}
	//}	
	//fclose(file_mesh_vtk);

	//cout << "Done writing mesh file file_mesh.vtk" << "\n";

	
	cout << "Done writing PLOT3D mesh file file_mesh.xyz" << "\n";

}


	// Create the boundary values for the domain, and send it to the function to create the mesh


void get_airfoil_coords(double*& x_af,double*& y_af,int n_afoil)
{
	// Airfoil coordinates (see http://airfoiltools.com/airfoil/naca4digit)
	double M = 0.0;
	double P = 0.0;
	double XX = 12.0;
	
	double m = M/100.0;
	double p = 0.1*P;
	double xx = XX/100.0;
	
	//double xUL[n_afoil], yUL[n_afoil];
	for (int j=0;j<n_afoil;++j)
	{	
		// Assign yc and dycdx
		double beta = pi - j * (2*pi)/(n_afoil-1);
		double xArf;
		if (abs(beta)>=pi/2.0){
			xArf = 1.0 - j*2.0/(n_afoil-1);					
		}
		else{
			xArf = (1.0 - cos(abs(beta)))/2.0;
		}

		xArf = abs(xArf);
		
		// Uncomment to cluster points at the leading edge
		//double xArf = (1.0 - cos(abs(beta)))/2.0;
				
		double yc, dycdx;
		if (xArf>=0 && xArf<p) {
			yc = m/pow(p,2.0) * (2.0*p*xArf - pow(xArf,2.0));
			dycdx = 2.0*m/pow(p,2.0) * (p - xArf);
		}
		else if (xArf>=p && xArf<=1.0){
			yc = m/pow(1.0-p,2.0) * (1.0 - 2.0*p + 2.0*p*xArf - pow(xArf,2.0));
			dycdx = 2.0*m/pow(1.0-p,2.0) * (p - xArf);
		}
		else {
			cout << "ERROR: xArf should not be greater than 1 or smaller than 0!" << endl;
			exit (EXIT_FAILURE);
		}
		
		double a0 = 0.2969;
		double a1 = -0.126;
		double a2 = -0.3516;
		double a3 = 0.2843;
		double a4 = -0.1036; //-0.1015;
		double yt = xx/0.2 * (a0*pow(xArf,0.5) + a1*xArf + a2*pow(xArf,2.0) + a3*pow(xArf,3.0) + a4*pow(xArf,4.0));
					
		double th = atan2(dycdx,1.0);
		
		if (beta>=0) {
			x_af[j] = xArf - yt*sin(th);
			y_af[j] = yc + yt*cos(th);					
		}
		else {
			x_af[j] = xArf + yt*sin(th);
			y_af[j] = yc - yt*cos(th);					
		}
	}
	
	// Write the airfoil coorfinates
	char filename[20];
	sprintf(filename, "airfoil_coords.dat");
	FILE *fid_err = fopen (filename, "w");
		
	for (int j=0;j<n_afoil;++j) {	
		fprintf(fid_err, "%e   %e\n",x_af[j],y_af[j]);
	}
	
	fclose(fid_err);
	cout << "===============================" << endl;
	printf("Wrote the airfoil coordinates in file: \"%s\" \n", filename);
	cout << "===============================" << endl;
	
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

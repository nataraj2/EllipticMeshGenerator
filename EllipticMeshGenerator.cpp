#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
using namespace std;

const double pi = M_PI;

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
void write_plot3d_grid(double***,double***,double***,int,int,int);
	
int main()
{
	int npts_afoil = 50;
	double* x_afoil	= new double[npts_afoil];
	double* y_afoil	= new double[npts_afoil];
	get_airfoil_coords(x_afoil,y_afoil,npts_afoil);
		
	int nx = 50;
	int ny = 25;
	int nz = 25;

	// Parameters defining the ellipsoid
	// (x-x0)^2/a^2 + (y-y0)^2/b^2 + (z-z_0)^2/c^2;

	double xmin = 0.0, xmax = 3.0, ymin = 0.0, ymax = 1.0, zmin = 0.0, zmax = 1.0;
	double a = 0.1, b = 0.3, c = 0.1;
	double x0 = 1.0, z0 = 0.5;

	// Number of smoothing iterations for the elliptic smoothing
	int n_iterations = 20;

	FILE* file_mesh_vtk;
	file_mesh_vtk = fopen("file_mesh.vtk", "w");
	fprintf(file_mesh_vtk,"%s\n", "# vtk DataFile Version 2.0");
	fprintf(file_mesh_vtk,"%s\n", "Simple example");
	fprintf(file_mesh_vtk,"%s\n", "ASCII");
	fprintf(file_mesh_vtk,"%s\n", "DATASET STRUCTURED_GRID");
	fprintf(file_mesh_vtk,"%s %d %d %d\n", "DIMENSIONS", nx, ny, nz);
	fprintf(file_mesh_vtk,"%s %d %s\n", "POINTS", nx*ny*nz, "float");

	// Define the coordinate arrays

	double*** x = Create3DMatrix(nx,ny,nz);
	double*** y = Create3DMatrix(nx,ny,nz);
	double*** z = Create3DMatrix(nx,ny,nz);

	//double x[nx][ny][nz], y[nx][ny][nz], z[nx][ny][nz];

	double delx = (xmax-xmin)/float(nx-1);	
	double dely = (ymax-ymin)/float(ny-1);	
	double delz = (zmax-zmin)/float(nz-1);

	cout << "The mesh dimensions is x, y, z are " << nx << " " << ny << " " << nz << "\n";

	cout << "Initializing mesh generation" << "\n";

	// Create the initial mesh slice-by-slice in the z direction
	for(int k=0;k<nz;k++){
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
				x[i][j][k] = 0.0;
				y[i][j][k] = 0.0;
				z[i][j][k] = 0.0;
			}
		}

		// The equation for the ellpsoid is ((x-x0)/a)^2 + ((y-y0)/b)^2 + ((z-z0)/c)^2 = 1

		// Stretch in z
		double delta = 4.0;
		double actual_pos = 0.5;
		double pos = 1.0-actual_pos;

		double zeta = float(k)/float(nz-1);
		double u1 = tanh(delta*(1.0-zeta))*(1.0-pos);
		double u2 = (2.0-tanh(delta*zeta))*pos;
		double fac = 1.0 - ((u1+u2)-pos);

		for(int i=0;i<nx; i++){
			for(int j=0;j<ny;j++){
				z[i][j][k] = zmin + fac*(zmax-zmin);
			}
		}	
	
		// Stretch in x
		double delta_x = 4.0;
		double actual_pos_x = 1-0.4;
		double pos_x = 1.0-actual_pos_x;

		double zval = z[0][0][k];

		// Assign values to the boundary values of x and y in this z slice
		if(fabs((zval-z0)/c)<1.0){
			for(int i=0;i<nx; i++){
	
				double xi = float(i)/float(nx-1);
				double u1 = tanh(delta_x*(1.0-xi))*(1.0-pos_x);
				double u2 = (2.0-tanh(delta_x*xi))*pos_x;
				double fac = 1.0 - ((u1+u2)-pos_x);

				double xval = xmin + fac*(xmax-xmin);
				x[i][0][k] = xval;
				x[i][ny-1][k] = xval;
				if(fabs((xval-x0)/a)<sqrt(1.0-(zval-z0)/c*(zval-z0)/c)){
					y[i][0][k] = b*sqrt(1.0-(xval-x0)/a*(xval-x0)/a - (zval-z0)/c*(zval-z0)/c);
				}
				else{
					y[i][0][k] = ymin;
				}
				y[i][ny-1][k] = ymax;
			}
		}
		else {
			for(int i=0;i<nx; i++){
				double xi = float(i)/float(nx-1);
				double u1 = tanh(delta_x*(1.0-xi))*(1.0-pos_x);
				double u2 = (2.0-tanh(delta_x*xi))*pos_x;
				double fac = 1.0 - ((u1+u2)-pos_x);

				double xval = xmin + fac*(xmax-xmin);

				x[i][0][k] = xval;
				x[i][ny-1][k] = xval;
				y[i][0][k] = ymin;
				y[i][ny-1][k] = ymax;
			}
		}


		// Stretch in y
		double delta_y = 2.0;

		for(int j=0;j<ny; j++){
			double eta = float(j)/float(ny-1);
			double fac = 1.0 + tanh(delta_y*(eta-1))/tanh(delta_y);
			double val = ymin + fac*(ymax-ymin);
			y[0][j][k] = val;
			y[nx-1][j][k] = val;
			x[0][j][k] = xmin;
			x[nx-1][j][k] = xmax;
		}

		// Do transfinite interpolation to fill the mesh in the domain interior using the boundary values of x and y defined above
		int m = nx-1;
		int n = ny-1;

		for(int i=1;i<nx-1; i++){
			for(int j=1;j<ny-1;j++){
				x[i][j][k] = float(i)/float(m)*x[m][j][k] + float(m-i)/float(m)*x[0][j][k] + float(j)/float(n)*x[i][n][k] + float(n-j)/float(n)*x[i][0][k] - float(i)/float(m)*float(j)/float(n)*x[m][n][k] - 
					  float(i)/float(m)*float(n-j)/float(n)*x[m][0][k] - float(m-i)/float(m)*float(j)/float(n)*x[0][n][k] - float(m-i)/float(m)*float(n-j)/float(n)*x[0][0][k];

			//y[i][j][k] = float(i)/float(m)*y[m][j][k] + float(m-i)/float(m)*y[0][j][k] + float(j)/float(n)*y[i][n][k] + float(n-j)/float(n)*y[i][0][k] - float(i)/float(m)*float(j)/float(n)*y[m][n][k] - 
			//		  float(i)/float(m)*float(n-j)/float(n)*y[m][0][k] - float(m-i)/float(m)*float(j)/float(n)*y[0][n][k] - float(m-i)/float(m)*float(n-j)/float(n)*y[0][0][k];


				double delta_y = 2.0;
				double eta = float(j)/float(ny-1);
				double fac = 1.0 + tanh(delta_y*(eta-1))/tanh(delta_y);
			
				y[i][j][k] = y[i][0][k] + fac*(y[i][n][k]-y[i][0][k]);
			}
		}
	
	} // End of k loop. Initial mesh generation for the whole domain ends here for the entire 3d domain

	cout << "Done with initial mesh generation" << "\n";


	// Do elliptic smoothing. Smooth only y. Elliptic equation is solved using Gauss-Siedel iteration
	for(int iter=0;iter<n_iterations;iter++){ 

		cout << "Doing smoothing iteration " << iter+1 << "\n";
		for(int i=1;i<nx-1;i++){
			for(int j=1;j<ny-1;j++){
				for(int k=1;k<nz-1;k++){

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

					a_11 = dxdxi*dxdxi     + dydxi*dydxi     + dzdxi*dzdxi;
					a_12 = dxdxi*dxdeta   + dydxi*dydeta    + dzdxi*dzdeta;
					a_13 = dxdxi*dxdzeta   + dydxi*dydzeta   + dzdxi*dzdzeta;
					a_22 = dxdeta*dxdeta   + dydeta*dydeta   + dzdeta*dzdeta;
					a_23 = dxdeta*dxdzeta  + dydeta*dydzeta  + dzdeta*dzdzeta;
					a_33 = dxdzeta*dxdzeta + dydzeta*dydzeta + dzdzeta*dzdzeta;

					alpha_11 = a_22*a_33 - a_23*a_23;
					alpha_12 = a_13*a_23 - a_12*a_33;
					alpha_13 = a_12*a_23 - a_13*a_22;
					alpha_22 = a_11*a_33 - a_13*a_13;
					alpha_23 = a_13*a_12 - a_11*a_23;
					alpha_33 = a_11*a_22 - a_12*a_12;
	
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
					
	
					double fac = 2.0*(alpha_11+alpha_22+alpha_33);
					
		
					//x[i][j] = (alpha*d2xdxi2_part - 2.0*beta*d2xdxideta + gamma*d2xdeta2_part)/fac;
					y[i][j][k] = (alpha_11*d2ydxi2_part + 2.0*alpha_12*d2ydxideta + alpha_22*d2ydeta2_part + alpha_33*d2ydzeta2_part )/fac;

				}
			}
		}
	}
			
	// Write VTK file
	for(int k=0; k<nz; k++){
       	for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				fprintf(file_mesh_vtk,"%g %g %g\n", x[i][j][k], y[i][j][k], z[i][j][k]);
			}
		}
	}	
	fclose(file_mesh_vtk);

	cout << "Done writing mesh file file_mesh.vtk" << "\n";

	// Write PLOT3D file
	write_plot3d_grid(x,y,z,nx,ny,nz);
	
	cout << "Done writing PLOT3D mesh file file_mesh.xyz" << "\n";


	// Write bc.dat

	// Find the ilimits of the geometry

	int istart, iend;
	for(int i=0; i<nx; i++){
		int j=0, k=0;
		if(x[i][j][k] < x0-a and x[i+1][j][k]>x0-a){
            istart = i;
        }
		if(x[i][j][k] < x0+a and x[i+1][j][k]>x0+a){
            iend = i+1;
        }
	}

	// Loop through the ilimits and then find the klimits for each i index (y is 0)

	FILE* file_bc_wall, *file_bc_farfield;
    file_bc_wall = fopen("file_bc_wall.vtk", "w");
    fprintf(file_bc_wall,"%s\n", "# vtk DataFile Version 2.0");
    fprintf(file_bc_wall, "%s\n", "Simple example");
    fprintf(file_bc_wall,"%s\n", "ASCII");
    fprintf(file_bc_wall,"%s\n", "DATASET POLYDATA");
	istart = istart+1;
	iend = iend-1;
    fprintf(file_bc_wall,"%s %d %s\n", "POINTS", (iend-istart+1)*2, "float");	


	file_bc_farfield = fopen("file_bc_farfield.vtk", "w");
    fprintf(file_bc_farfield,"%s\n", "# vtk DataFile Version 2.0");
    fprintf(file_bc_farfield, "%s\n", "Simple example");
    fprintf(file_bc_farfield,"%s\n", "ASCII");
    fprintf(file_bc_farfield,"%s\n", "DATASET POLYDATA");
    fprintf(file_bc_farfield,"%s %d %s\n", "POINTS", (iend-istart+1)*4, "float");	


	FILE *bc_file;
	bc_file = fopen("bc.dat","w");
	fprintf(bc_file, "%s\n", "#Bottom Wall");

	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 24, 2, 1, istart-1, 1, 1, 1, -1);
	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 24, 2, iend+1, -1, 1, 1, 1, -1);

	for(int i=istart; i<=iend; i++){
		// Find the z limits (x-x0)^2/a^2 + (z-z0)^2/c^2 = 1

		double xval = x[i][0][0];
		
		double len = c*sqrt((1.0 - (xval-x0)*(xval-x0)/(a*a)));
		double zmin = z0 - len;
		double zmax = z0 + len;
		
		int kstart, kend;
		for(int k=0; k<nz; k++){
        	int j = 0;
        	if(z[i][j][k] < z0-len and z[i][j][k+1]>z0-len){
            	kstart = k;
        	}
        	if(z[i][j][k] < z0+len and z[i][j][k+1]>z0+len){
            	kend = k+1;
        	}
		}
		// Write into wall vtk
		fprintf(file_bc_wall,"%g %g %g\n", x[i][0][0], 0.0, z[0][0][kstart]);	
		fprintf(file_bc_wall,"%g %g %g\n", x[i][0][0], 0.0, z[0][0][kend]);

		// Write into far-field vtk
		fprintf(file_bc_farfield,"%g %g %g\n", x[i][0][0], 0.0, z[0][0][0]);	
		fprintf(file_bc_farfield,"%g %g %g\n", x[i][0][0], 0.0, z[0][0][kstart]);	
		fprintf(file_bc_farfield,"%g %g %g\n", x[i][0][0], 0.0, z[0][0][kend]);	
		fprintf(file_bc_farfield,"%g %g %g\n", x[i][0][0], 0.0, z[0][0][nz-1]);	

		// Write into bc.dat
		fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 22, 2, i, i, 1, 1, kstart, kend);
		fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 24, 2, i, i, 1, 1, 1, kstart);
		fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 24, 2, i, i, 1, 1, kend, -1);
		
	}

		fprintf(file_bc_wall,"%s %d %d\n", "LINES", (iend-istart+1), (iend-istart+1)*3);
		for(int num=0;num<(iend-istart+1);num++){
    		fprintf(file_bc_wall,"%d %d %d\n", 2, 2*num, 2*num+1);
		}
		fprintf(file_bc_farfield,"%s %d %d\n", "LINES", (iend-istart+1)*2, (iend-istart+1)*2*3);
		for(int num=0;num<(iend-istart+1);num++){
    		fprintf(file_bc_farfield,"%d %d %d\n", 2, 4*num, 4*num+1);
    		fprintf(file_bc_farfield,"%d %d %d\n", 2, 4*num+2, 4*num+3);
		}

	fclose(file_bc_wall);
	fclose(file_bc_farfield);
	fclose(bc_file);

	cout << "Done writing boundary conditions file bc.dat" << "\n";
	// Find the limits

	// Find the i index of the center of the hump

	/*int istart, iend, kstart, kend;
	for(int i=0; i<nx; i++){
	int j = 0, k = 0;	
        if(x[i][j][k] < x0-a and x[i+1][j][k]>x0-a){
			istart = i;
            cout << "istart is " << i << "\n";
        }
		if(x[i][j][k] < x0+a and x[i+1][j][k]>x0+a){
			iend = i+1;
            cout << "iend is " << i << "\n";
        }
    }

		for(int k=0; k<nz; k++){
		int i = 0, j = 0;
        if(z[i][j][k] < z0-c and z[i][j][k+1]>z0-c){
			kstart = k;
            cout << "kstart is " << k << "\n";
        }
		if(z[i][j][k] < z0+c and z[i][j][k+1]>z0+c){
			kend = k+1;
            cout << "kend is " << k << "\n";
        }
    }


	// Write vtk file for lines of the boundary

	FILE* file_bc_lines;
	file_bc_lines = fopen("file_bc_lines.vtk", "w");
	fprintf(file_bc_lines,"%s\n", "# vtk DataFile Version 2.0");
	fprintf(file_bc_lines,"%s\n", "Simple example");
	fprintf(file_bc_lines,"%s\n", "ASCII");
	fprintf(file_bc_lines,"%s\n", "DATASET POLYDATA");
	fprintf(file_bc_lines,"%s %d %s\n", "POINTS", 4, "float");
	fprintf(file_bc_lines,"%g %g %g\n", x[istart][0][0], 0.0, z[0][0][kstart]);	
	fprintf(file_bc_lines,"%g %g %g\n", x[iend][0][0], 0.0, z[0][0][kstart]);	
	fprintf(file_bc_lines,"%g %g %g\n", x[iend][0][0], 0.0, z[0][0][kend]);	
	fprintf(file_bc_lines,"%g %g %g\n", x[istart][0][0], 0.0, z[0][0][kend]);

	fprintf(file_bc_lines,"%g %g %g\n", x[istart][0][0], 0.0, z[0][0][0]);
	fprintf(file_bc_lines,"%g %g %g\n", x[iend][0][0], 0.0, z[0][0][0]);
	fprintf(file_bc_lines,"%g %g %g\n", x[istart][0][0], 0.0, z[0][0][nz-1]);
	fprintf(file_bc_lines,"%g %g %g\n", x[istart][0][0], 0.0, z[0][0][nz-1]);

	fprintf(file_bc_lines,"%s %d %d\n", "LINES", 4, 12);
	fprintf(file_bc_lines,"%d %d %d\n", 2, 0, 1);
	fprintf(file_bc_lines,"%d %d %d\n", 2, 1, 2);
	fprintf(file_bc_lines,"%d %d %d\n", 2, 2, 3);
	fprintf(file_bc_lines,"%d %d %d\n", 2, 3, 0);

	fclose(file_bc_lines);


	// Write bc.dat

	FILE *bc_file;
	bc_file = fopen("bc.dat","w");
	fprintf(bc_file, "%s\n", "#Bottom Wall");
	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 22, 2, istart, iend, 1, 1, kstart, kend);
	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 21, 2, istart, iend, 1, 1, 1, kstart);
	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 21, 2, istart, iend, 1, 1, kend, -1);
	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 21, 2, 1, istart, 1, 1, 1, -1);
	fprintf(bc_file, "%d %d %d %d %d %d %d %d %d\n", 1, 21, 2, iend, -1, 1, 1, 1, -1);

	fclose(bc_file);*/
	
	Delete3DMatrix(x, nx, ny, nz);
	Delete3DMatrix(y, nx, ny, nz);
	Delete3DMatrix(z, nx, ny, nz);
	
	return 0;

}

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

void write_plot3d_grid(double*** x,double*** y,double*** z,int Nx,int Ny,int Nz)
{
	ofstream outFile;
	char filename[] = "test_grid3D.xyz";
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

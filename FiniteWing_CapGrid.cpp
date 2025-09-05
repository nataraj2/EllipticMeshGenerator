#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

static constexpr double PI = 3.14159265358979323846 ;
static constexpr double pi = 3.14159265358979323846 ;
struct Grid {
    int Nx, Ny, Nz;
    double*** x;
    double*** y;
    double*** z;
};
inline double clamp01(double x){ return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }

inline double naca0012_thickness(double x){
    x = std::max(0.0, std::min(1.0, x));
    const double t = 0.12; 
    const double a0 =  0.2969;
    const double a1 = -0.1260;
    const double a2 = -0.3516;
    const double a3 =  0.2843;
    const double a4 = -0.1015; 
    const double sqrtx = std::sqrt(std::max(0.0, x));
    const double yt = 5.0 * t * (a0*sqrtx + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x);
    return std::fabs(yt);
}
double*** Create3DMatrix(int nx, int ny, int nz) {
    double*** matrix = new double**[nx];
    for(int i = 0; i < nx; i++) {
        matrix[i] = new double*[ny];
        for(int j = 0; j < ny; j++) {
            matrix[i][j] = new double[nz];
        }
    }
    return matrix;
}
int*** Create3DMatrix_INT(int nx, int ny, int nz)
{
    int*** matrix = new int**[nx];
    for(int i = 0; i < nx; i++) {
        matrix[i] = new int*[ny];
        for(int j = 0; j < ny; j++) {
            matrix[i][j] = new int[nz];
        }
    }
    return matrix;
}
void Delete3DMatrix(double*** matrix, int nx, int ny, int nz) {
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            delete[] matrix[i][j];
        }
        delete[] matrix[i];
    }
    delete[] matrix;
}

void create_mesh(int nx, int ny, int nz, double***& x, double***&y, double***& z, double zmin, double zmax, 
				 const double *x_xi_0, const double *x_xi_max, const double *x_eta_0, const double *x_eta_max,
				 const double *y_xi_0, const double *y_xi_max, const double *y_eta_0, const double *y_eta_max, double delta_y, 
				 int n_iterations, bool invert, std::string filename)
{
	std::cout << "The mesh dimensions is x, y, z are " << nx << " " << ny << " " << nz << "\n";

	std::cout << "Initializing mesh generation" << "\n";	

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
		double delz;
		if (nz>1) {
			double zeta = float(k)/float(nz-1);
			double u1 = tanh(delta*(1.0-zeta))*(1.0-pos);
			double u2 = (2.0-tanh(delta*zeta))*pos;
			double fac = 1.0 - ((u1+u2)-pos);
		
			delz = (zmax-zmin)/(nz-1);
		}
		else {
			delz = 0.0;
		}

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

	std::cout << "Done with initial mesh generation" << "\n";


	// Do elliptic smoothing. Smooth only y. Elliptic equation is solved using Gauss-Siedel iteration
	for(int iter=0;iter<n_iterations;iter++){ 

		std::cout << "Doing smoothing iteration " << iter+1 << "\n";
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

	// write_plot3d_grid(x,y,z,nx,ny,nz,"mesh_OGrid.xyz");
	
	/*gl->Nx = nx;
	gl->Ny = ny;
	gl->Nz = nz;
	gl->x = Create3DMatrix(nx,ny,nz);
	gl->y = Create3DMatrix(nx,ny,nz);
	gl->z = Create3DMatrix(nx,ny,nz);*/
	for (int k=0;k<nz;++k)
		for (int j=0;j<ny;++j)
			for (int i=0;i<nx;++i)
			{
				gl-> x[j][(nx-1)-i][k]= x[i][j][k];
				gl-> y[j][(nx-1)-i][k]= y[i][j][k];
				gl->z[j][(nx-1)-i][k] = z[i][j][k];
			}
			
	Delete3DMatrix(x, nx, ny, nz);
	Delete3DMatrix(y, nx, ny, nz);
	Delete3DMatrix(z, nx, ny, nz);
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
			std::cout << "ERROR: xArf should not be greater than 1 or smaller than 0!" << std::endl;
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

	
}

int main(int argc, char** argv){

////cap grid generation/////// 
    int Nr = 20;        
    int Ntheta = 41;     
    int kSkip = 3; //added to generate the first slice exactly matching LE
    int Nz = 20 + kSkip;      //newNz in line 47 is the actual slice no  

    
    int arcSlices = 8; //no of arc slices to have
    int jClamp = 1; //to have orthogonality (1 j row)
    double r_blend = 0.01;
    double fracO = 0.9;      //fraction of Ntheta for Ogrid portion of those extended slices (extended = semi ellipse + two rectangles at bottom)
    double RotateAngle = 0.0; 

    const int newNz = Nz - kSkip;

    auto make3d = [&](int K, int J, int I){
        return vector<vector<vector<double>>>(K, vector<vector<double>>(J, vector<double>(I, 0.0)));
    };
    auto X = make3d(newNz, Ntheta, Nr);
    auto Y = make3d(newNz, Ntheta, Nr);
    auto Z = make3d(newNz, Ntheta, Nr);

    ///func to build those extended slices
    auto build_extended_slice = [&](double a_ell, double b_ell, double R_outer,
                                          vector<vector<double>>& Y2d,
                                          vector<vector<double>>& Z2d,
                                          vector<double>& tiltW,
                                          double z_value)
    {

        vector<vector<double>> Y1(Ntheta, vector<double>(Nr, 0.0));
        vector<vector<double>> Z1(Ntheta, vector<double>(Nr, 0.0));
        double rect_start = -0.3;

        for(int j=0;j<Ntheta;++j){
            double theta = PI * j / (Ntheta - 1);
            double cosT = cos(theta), sinT = sin(theta);
            double yin = a_ell * cosT;
            double zin = b_ell * sinT;
            double yout = R_outer * cosT;
            double zout = R_outer * sinT;
            for(int i=0;i<Nr;++i){
                double s = (Nr==1)?0.0: double(i)/(Nr-1);
                Y1[j][i] = (1.0 - s)*yin + s*yout;
                Z1[j][i] = (1.0 - s)*zin + s*zout;
            }
        }

        //Left rectangle
        vector<vector<double>> Y2(Nr, vector<double>(Nr, 0.0));
        vector<vector<double>> Z2(Nr, vector<double>(Nr, 0.0));
        for(int j=0;j<Nr;++j){
            double y = -R_outer + (R_outer - a_ell) * j / (Nr - 1);
            for(int i=0;i<Nr;++i){
                double z = ((b_ell - z_value) * 2.9 * i / (Nr - 1)) +rect_start; //2.9 is multiplied to increase rectangular dz spacing, more to increase rect spacing
                Y2[j][i] = y;
                Z2[j][i] = z + z_value;
            }
        }
        //Right rectangle 
        auto Y3 = Y2, Z3 = Z2;
        for(int j=0;j<Nr;++j){
            double y = a_ell + (R_outer - a_ell) * j / (Nr - 1);
            for(int i=0;i<Nr;++i){
                double z = ((b_ell - z_value) * 2.9 * i / (Nr - 1)) +rect_start; 
                Y3[j][i] = y;
                Z3[j][i] = z + z_value;
            }
        }

        //Merge along j
        const int total_j = Nr + (Ntheta - 1) + (Nr - 1);
        vector<vector<double>> Y_full(total_j, vector<double>(Nr, 0.0));
        vector<vector<double>> Z_full(total_j, vector<double>(Nr, 0.0));

        for(int j=0;j<Nr;++j) for(int i=0;i<Nr;++i){ Y_full[j][i] = Y3[i][j]; Z_full[j][i] = Z3[i][j]; }
        //O-grid
        for(int j=1;j<Ntheta;++j){ int jg = Nr - 1 + j; for(int i=0;i<Nr;++i){ Y_full[jg][i] = Y1[j][i]; Z_full[jg][i] = Z1[j][i]; }}
        //Left rectangle (flip; skip j=0)
        for(int j=1;j<Nr;++j){ int jg = (Nr - 1) + (Ntheta - 1) + j; for(int i=0;i<Nr;++i){ int iF = Nr-1-i, jF = Nr-1-j; Y_full[jg][i] = Y2[iF][jF]; Z_full[jg][i] = Z2[iF][jF]; }}

        //Resample to Ntheta with symmetric rectangle distribution
        Y2d.assign(Ntheta, vector<double>(Nr, 0.0));
        Z2d.assign(Ntheta, vector<double>(Nr, 0.0));
        int Co = std::max(1, std::min(Ntheta - 2, (int)std::lround(fracO * (Ntheta - 2))));
        int rem = Ntheta - Co;
        if(rem < 2){ Co = std::max(1, Ntheta - 2); rem = Ntheta - Co; }
        if(rem % 2 != 0){ if(Co < Ntheta - 2) ++Co; else --Co; rem = Ntheta - Co; }
        int Cr = rem/2, Cl = rem/2;
        const int srcR_start = 0,            srcR_count = std::max(1, Nr-1),       srcR_last = srcR_start + srcR_count - 1;
        const int srcO_start = (Nr - 1) + 1, srcO_count = std::max(1, Ntheta-2),    srcO_last = srcO_start + srcO_count - 1;
        const int srcL_start = (Nr - 1) + (Ntheta - 1) + 1,
                  srcL_count = std::max(1, Nr-1),
                  srcL_last  = srcL_start + srcL_count - 1;

        auto sample_block = [&](int C, int s0, int s1, int jlocal, int i){
            int L = std::max(1, s1 - s0);
            double s = (C==1)? 0.0 : double(jlocal)/double(C-1);
            double jf = double(s0) + s * double(L);
            int j0 = (int)std::floor(jf);
            int j1 = std::min(j0 + 1, s1);
            double t = jf - j0;
            double y = (1.0 - t) * Y_full[j0][i] + t * Y_full[j1][i];
            double z = (1.0 - t) * Z_full[j0][i] + t * Z_full[j1][i];
            return std::pair<double,double>(y,z);
        };

        for(int j=0;j<Ntheta;++j){
            for(int i=0;i<Nr;++i){
                if(j < Cr){ auto yz = sample_block(Cr, srcR_start, srcR_last, j, i); Y2d[j][i]=yz.first; Z2d[j][i]=yz.second; }
                else if(j < Cr + Co){ int jl = j - Cr; auto yz = sample_block(Co, srcO_start, srcO_last, jl, i); Y2d[j][i]=yz.first; Z2d[j][i]=yz.second; }
                else { int jl = j - (Cr + Co); auto yz = sample_block(Cl, srcL_start, srcL_last, jl, i); Y2d[j][i]=yz.first; Z2d[j][i]=yz.second; }
            }
        }
        //first and last j to be orthogonal, weight
        tiltW.assign(Ntheta, 0.0);
        auto smoothstep = [&](double x){ x = clamp01(x); return x*x*(3.0 - 2.0*x); };
        const double m_blend = 3.0;
        for(int j=0;j<Ntheta;++j){ double w = smoothstep(double(j)/m_blend); if(j==0 || j==Ntheta-1) w=0.0; tiltW[j] = std::min(1.0, std::max(0.0, w)); }
    };

    // LE/TE z profiles to follow the arc, needs cleaning
    static bool le_params_init = false;
    vector<double> z_le_by_k, z_te_by_k, bell_by_k;
    if(arcSlices > 0){
        z_le_by_k.resize(arcSlices);
        z_te_by_k.resize(arcSlices);
        bell_by_k.assign(arcSlices, 0.15);
        const double z0 = 0.03;  //baseline
        const double dz = 0.08;  //total rise
        auto smoother = [&](double t){ t = clamp01(t); return t*t*t*(t*(t*6.0 - 15.0) + 10.0); };
        for(int kk=0; kk<arcSlices; ++kk){
            double t = (arcSlices>1)? double(kk)/double(arcSlices-1) : 0.0;
            double u = smoother(t);
            z_le_by_k[kk] = z0 + dz * u;           //rise
            z_te_by_k[kk] = z0 + dz * (1.0 - u);   //mirror fall
        }
        if(arcSlices > 2) bell_by_k[2] = 0.105;
        if(arcSlices > 3) bell_by_k[3] = 0.125;
        le_params_init = true;
    }

    int outK = 0;  //k slices will loop over outK, basically outK=0 is first slices (LE), outK= newNz -1  (or newNzLoc -1 )is last (TE)
    for(int k=0; k<Nz; ++k){
        if(k < kSkip) continue; //skip a few at start to fit LE slice
        double spacing  = double(k) / double(Nz - 1);
        double x_c   = spacing; 
        double thick = naca0012_thickness(x_c);
        double r_in  = std::max(thick, r_blend);
        double r_out = 0.12 * 8.5;
        double dtheta = PI / (Ntheta - 1);

        //slice angles & conditions(lots of H-codings)
        double phi_lead = 0.0, phi_trail = 0.0;
        const int newNzLoc = newNz; //new
        if(outK < arcSlices && arcSlices>1){
            phi_lead = -0.5*PI + (double(outK)/double(arcSlices-1))*(0.5*PI);
        }
        if(outK >= newNzLoc - arcSlices && arcSlices>1){
            int tIdx = outK - (newNzLoc - arcSlices);
            phi_trail = (double(tIdx)/double(arcSlices-1))*(0.5*PI);
        }

        bool isRigidFirst  = (outK == 0);
        bool isRigidLast   = (outK == newNzLoc - 1);
        bool isRigidSecond = (outK == 1);
        bool isRigidSecondLast = (outK == newNzLoc - 2);
        bool isRigidThirdLast  = (outK == newNzLoc - 3);
        if(isRigidSecond)      r_out = 0.12 * 9.0;
        if(isRigidSecondLast)  r_out = 0.12 * 10.0;
        if(isRigidLast)        r_out = 0.12 * 10.0;
        bool use_extended_slice = (!isRigidFirst && !isRigidLast && !isRigidSecond && !isRigidSecondLast);

        double z_value = 0.0;
        double b_ell_cur = 0.15;
        if(arcSlices > 0){
            if(outK < arcSlices){ z_value = z_le_by_k[outK]; b_ell_cur = bell_by_k[outK]; }
            else if(outK >= newNzLoc - arcSlices){ int tIdx = outK - (newNzLoc - arcSlices); z_value = z_te_by_k[tIdx]; }
            else { z_value = z_le_by_k.back(); b_ell_cur = 0.15; }
        }

        vector<vector<double>> Ycode2, Zcode2; vector<double> tiltW;
        if(use_extended_slice){
            const double tNACA = std::max(naca0012_thickness(x_c), r_blend);
            double a_ell = tNACA;
            double b_ell = 0.15;
            double R_outer = 0.12 * 10.0;
            double z_val_local = 0.0;
            if(arcSlices > 0){
                if(outK < arcSlices){ z_val_local = z_le_by_k[outK]; b_ell = bell_by_k[outK]; }
                else if(outK >= newNzLoc - arcSlices){ int tIdx = outK - (newNzLoc - arcSlices); z_val_local = z_te_by_k[tIdx]; }
                else { z_val_local = z_le_by_k.back(); }
            }
            if(outK == 2){ a_ell = a_ell - 0.016; R_outer = 0.12 * 9.5; }
            if(outK == newNzLoc - 3){ b_ell = b_ell - 0.03; }
            if(outK == newNzLoc - 4){ b_ell = b_ell + 0.001; }
            if(outK == newNzLoc - 5){ b_ell = b_ell + 0.001; }
            build_extended_slice(a_ell, b_ell, R_outer, Ycode2, Zcode2, tiltW, z_val_local);
        }

        //per-j loop to compute coordinates
        for(int j=0; j<Ntheta; ++j){
            double theta = j * dtheta;
            bool atClamp = (j <= jClamp) || (j >= Ntheta - 1 - jClamp);
            double local_tilt = 0.0;
            if(isRigidFirst){ local_tilt = -0.5 * PI; }
            else if(isRigidSecond){ local_tilt = -0.41 * PI; }
            else if(outK == 2){ local_tilt = -0.3 * PI; }
            else if(outK == 3){ local_tilt = -0.2 * PI * std::sin(theta); }
            else if(outK < arcSlices && !atClamp){ local_tilt = phi_lead * std::sin(theta); }
            else if(isRigidLast){ local_tilt = 0.5 * PI; }
            else if(isRigidSecondLast){ local_tilt = 0.41 * PI; }
            else if(isRigidThirdLast){ local_tilt = 0.35 * PI; }
            else if(outK >= newNzLoc - arcSlices && !atClamp){ local_tilt = phi_trail * std::sin(theta); }

            double tilt_mask = 1.0;
            if(use_extended_slice) tilt_mask = tiltW[j];
            if(atClamp) tilt_mask = 0.0;
            double local_tilt_effective = local_tilt * tilt_mask;

            if(isRigidFirst || isRigidSecond){
                //LE shape profile
                double RX_ANGLE = -PI/2.0; const double c_rx = cos(RX_ANGLE), s_rx = sin(RX_ANGLE);
                double R_outer = 0.12 * 9.2; const double x_max = 0.04;
                for(int jj=0; jj<Ntheta; ++jj){
                    const double th = PI * double(jj) / double(Ntheta - 1);
                    double xin, yin;
                    if(th <= 0.5 * PI){ double s01 = th / (0.5 * PI); xin = x_max * (1.0 - s01); yin = +naca0012_thickness(xin); }
                    else { double s01 = (th - 0.5 * PI) / (0.5 * PI); xin = x_max * s01; yin = -naca0012_thickness(xin); }
                    const double xout_xy = -R_outer * sin(th);
                    const double yout_xy =  R_outer * cos(th);
                    const double Yin = xin;
                    const double Zin = yin;
                    double Yout = xout_xy; const double Zout = yout_xy;
                    if(jj==0 || jj==Ntheta-1) Yout = Yin;
                    for(int i=0;i<Nr;++i){
                        const double s = (Nr>1)? double(i)/double(Nr-1) : 0.0;
                        const double y0 = (1.0 - s) * Yin + s * Yout;
                        const double z0 = (1.0 - s) * Zin + s * Zout;
                        const double Yx = c_rx * y0 - s_rx * z0;
                        const double Zx = s_rx * y0 + c_rx * z0;
                        const double st = sin(local_tilt);
                        const double ct = cos(local_tilt);
                        double Xv = x_c + Zx * st;
                        double Yv = Yx;
                        double Zv = Zx * ct;
                        if(isRigidFirst){ Xv += -0.0938 - 0.04 - 0.002; }
                        if(isRigidSecond && outK>0){ Xv = X[outK-1][jj][i]; Zv += 0.08; }
                        X[outK][jj][i] = Xv;
                        Y[outK][jj][i] = Yv;
                        Z[outK][jj][i] = Zv - 0.3;
                    }
                }
            }
            else if(use_extended_slice){
                //inner radious i=0 with amplified tilt to have a smoother arc, needs to be fixed
                {
                    const int i0 = 0;
                    double y0 = Ycode2[j][i0];
                    double z0 = Zcode2[j][i0];
                    double dz01 = (Nr>1)? (Zcode2[j][1] - Zcode2[j][0]) : 1e-9;
                    double alpha = 1.5;
                    double z_pivot = z0 - alpha * dz01;
                    double TILT_GAIN = 1.35;
                    double st = sin(TILT_GAIN * local_tilt_effective);
                    double ct = cos(TILT_GAIN * local_tilt_effective);
                    double Xv0 = x_c + (z0 - z_pivot) * st;
                    double Yv0 = y0;
                    double Zv0 = z_pivot + (z0 - z_pivot) * ct;
                    X[outK][j][i0] = Xv0; Y[outK][j][i0] = Yv0; Z[outK][j][i0] = Zv0;
                }
                for(int i=1;i<Nr;++i){
                    double y0 = Ycode2[j][i];
                    double z0 = Zcode2[j][i];
                    double dz01 = (Nr>1)? (Zcode2[j][1] - Zcode2[j][0]) : 1e-9;
                    double alpha = 1.5;
                    double z_pivot = Zcode2[j][0] - alpha * dz01;
                    double TILT_GAIN = 1.0;
                    double st = sin(TILT_GAIN * local_tilt_effective);
                    double ct = cos(TILT_GAIN * local_tilt_effective);
                    double Xv = x_c + (z0 - z_pivot) * st;
                    double Yv = y0;
                    double Zv = z_pivot + (z0 - z_pivot) * ct;
                    X[outK][j][i] = Xv; Y[outK][j][i] = Yv; Z[outK][j][i] = Zv;
                }
            }
            else if(isRigidSecondLast){
                for(int i=0;i<Nr;++i){
                    double r_in_theta = std::max(r_in, r_blend);
                    double dr_local = (r_out - r_in_theta) / std::max(1, Nr - 1);
                    double r = r_in_theta + i * dr_local;
                    double y0 = r * cos(theta);
                    double z0 = r * sin(theta);
                    double Xv = x_c + z0 * sin(local_tilt);
                    double Yv = y0;
                    double Zv = z0 * cos(local_tilt);
                    Zv += 0.08 - 0.011;
                    X[outK][j][i] = Xv; Y[outK][j][i] = Yv; Z[outK][j][i] = Zv - 0.3;
                }
            }
            else if(isRigidLast){
                const double r_outer = 0.12 * 10.0;
                const double tip_x = 0.05, tip_y = 0.0;
                const double base_low_x = 0.0, base_low_y = -0.01 + 0.0028;
                const double base_up_x  = 0.0, base_up_y  =  0.01 - 0.0028;
                vector<double> xb(Ntheta), yb(Ntheta), xo(Ntheta), yo(Ntheta);
                const double dth = PI / (Ntheta - 1);
                const int j_mid = (Ntheta - 1)/2;
                for(int jj=0; jj<Ntheta; ++jj){
                    double theta2 = PI/2 - jj * dth;
                    if(jj <= j_mid){ double uL = (j_mid>0)? double(jj)/double(j_mid) : 0.0; xb[jj] = base_up_x + uL*(tip_x - base_up_x);  yb[jj] = base_up_y + uL*(tip_y - base_up_y); }
                    else { double denom = double((Ntheta-1)-j_mid); double uU = (denom>0.0)? double(jj - j_mid)/denom : 1.0; xb[jj] = tip_x + uU*(base_low_x - tip_x); yb[jj] = tip_y + uU*(base_low_y - tip_y); }
                    xo[jj] = r_outer * cos(theta2);
                    yo[jj] = r_outer * sin(theta2);
                }
                for(int jj=0; jj<Ntheta; ++jj){
                    double local_tilt2 = 0.5*PI;
                    for(int i=0;i<Nr;++i){
                        double s = double(i)/double(Nr-1);
                        double y0 = (1 - s) * yb[jj] + s * yo[jj];
                        double z0 = (1 - s) * xb[jj] + s * xo[jj];
                        double Xv = x_c + z0 * sin(local_tilt2);
                        double Yv = y0;
                        double Zv = z0 * cos(local_tilt2);
                        X[outK][jj][i] = Xv; Y[outK][jj][i] = Yv; Z[outK][jj][i] = Zv - 0.3;
                    }
                }
            }
        } 
        //end j loop

        //function to shift,rotate a little bit if needed
        auto rotate_slice_about_innerY = [&](int kk, double rot_deg, double xShift, double zShift){
            double ang = rot_deg * PI / 180.0; double c = cos(ang), s = sin(ang);
            for(int j2=0;j2<Ntheta;++j2){
                double x0 = X[kk][j2][0];
                double z0 = Z[kk][j2][0];
                for(int i2=0;i2<Nr;++i2){
                    double Xp = X[kk][j2][i2] - x0;
                    double Zp = Z[kk][j2][i2] - z0;
                    double Xrot =  c*Xp + s*Zp;
                    double Zrot = -s*Xp + c*Zp;
                    X[kk][j2][i2] = Xrot + x0 + xShift;
                    Z[kk][j2][i2] = Zrot + z0 + zShift;
                }
            }
        };
        //H-coded, needs to be fixed
        if(outK == 2) rotate_slice_about_innerY(outK, -20.0, -0.16, 0.07);
        if(outK == 3) rotate_slice_about_innerY(outK, -20.0, -0.08, 0.02);
        if(outK == 4) rotate_slice_about_innerY(outK,  -0.0, -0.04, 0.02);
        if(outK == 5) rotate_slice_about_innerY(outK,  -0.0, -0.00, 0.01);
        if(outK == 6) rotate_slice_about_innerY(outK,  -0.0, -0.00, 0.005);
        if(outK == newNzLoc - 3) rotate_slice_about_innerY(outK, +10.0, 0.0, 0.04);
        if(outK == newNzLoc - 4) rotate_slice_about_innerY(outK, +10.0, 0.0, 0.02);
        if(outK == newNzLoc - 5) rotate_slice_about_innerY(outK,  +0.0, 0.0, 0.01);

        ++outK;
    } //end k loop (build slices)

    //H-coded, needs to be fixed
    {
        const int k_last  = newNz - 1;
        const int k_slast = newNz - 2;
        const int k_tlast = newNz - 3;
        for(int j=0;j<Ntheta;++j){
            for(int i=0;i<Nr;++i){
                X[k_last][j][i]  = X[k_last][j][i] - 0.05;
                X[k_slast][j][i] = X[k_last][j][i];
                X[k_tlast][j][i] = X[k_tlast][j][i] + 0.02;
                Z[k_tlast][j][i] = Z[k_tlast][j][i] + 0.04;
            }
        }
    }

    //for bottom arc 
    auto smootherstep5 = [&](double t){ t = clamp01(t); return t*t*t*(t*(t*6.0 - 15.0) + 10.0); };
    vector<double> x_rep(newNz, 0.0);
    for(int k=0;k<newNz;++k){ double sumx = 0.0; for(int j=0;j<Ntheta;++j) sumx += X[k][j][0]; x_rep[k] = sumx / std::max(1, Ntheta); }
    double xmin = x_rep[0], xmax = x_rep[0];
    for(double v: x_rep){ xmin = std::min(xmin, v); xmax = std::max(xmax, v); }
    double dx = xmax - xmin; bool ok = (dx > 1e-14);
    const double s_le_end = 0.12, s_te_beg = 0.88;
    const double z_base = 0.0, z_peak = 0.08;
    auto z_profile = [&](double s){ s = clamp01(s); if(s <= s_le_end){ double u = (s_le_end>0.0)? (s/s_le_end) : 1.0; return z_base + (z_peak - z_base) * smootherstep5(u); } if(s < s_te_beg){ return z_peak; } double u = (1.0 - s_te_beg > 0.0)? ((s - s_te_beg)/(1.0 - s_te_beg)) : 1.0; return z_peak + (z_base - z_peak) * smootherstep5(u); };
    for(int k=0;k<newNz;++k){ double s = ok? ((x_rep[k] - xmin)/dx) : (double(k)/std::max(1, newNz-1)); double z_off = z_profile(s); if(k==0 || k==newNz-1) z_off = z_base; for(int j=0;j<Ntheta;++j) for(int i=0;i<Nr;++i) Z[k][j][i] += z_off; }

    //filling coordinate array 
    {
        double y_angle = RotateAngle * PI / 180.0; // yaw about Y
        double z_angle = 0.0;                      // keep as in your code
        double cos_y = cos(y_angle), sin_y = sin(y_angle);
        double cos_z = cos(z_angle), sin_z = sin(z_angle);
        // Also the original code lifted Z by +1.98; we mirror that before rotation
        for(int k=0;k<newNz;++k){
            for(int j=0;j<Ntheta;++j){
                for(int i=0;i<Nr;++i){
                    double x_orig = X[k][j][i];
                    double y_orig = Y[k][j][i];
                    double z_orig = Z[k][j][i] + 1.98;
                    double x1 = x_orig * cos_z - y_orig * sin_z;
                    double y1 = x_orig * sin_z + y_orig * cos_z;
                    double z1 = z_orig;
                    double x2 = x1 * cos_y + z1 * sin_y;
                    double y2 = y1;
                    double z2 = -x1 * sin_y + z1 * cos_y;
                    X[k][j][i] = x2; 
                    Y[k][j][i] = y2; 
                    Z[k][j][i] = z2;
                }
            }
        }
    }

/////end of grid generation for cap grid/////
    
// wing grid generation

    Grid temp_grid;

    const int Ni = 30;
    const int Nj = 60;
    const int Nk = 20;

  //O-grid parameters
    double flat_portion = 15.0;        
    double semicircle_portion = 20.0;
    double x_circle_start = 0.0;       
    double len_y = 2.0;               
    double delta_y = 1.2;             
    int n_iterations = 20;           
    //Get airfoil coordinates
    int npts_afoil = Nj;
    double* x_afoil = new double[npts_afoil];
    double* y_afoil = new double[npts_afoil];
    get_airfoil_coords(x_afoil, y_afoil, npts_afoil);
    int indx[2] = {0, npts_afoil-1};

    temp_grid.Nx = Ni; 
    
    temp_grid.Ny = Nj; 
    temp_grid.Nz = Nk; 
    temp_grid.x = Create3DMatrix(temp_grid.Nx, temp_grid.Ny, temp_grid.Nz);
    temp_grid.y = Create3DMatrix(temp_grid.Nx, temp_grid.Ny, temp_grid.Nz);
    temp_grid.z = Create3DMatrix(temp_grid.Nx, temp_grid.Ny, temp_grid.Nz);
    double zs[] = {0, -1.5}; //if i want to rotate with respect to 0 point, then set 0 here
	double ze[] = {2, 1.5};
    //Generate O-grid
    create_OGrid_mesh(&temp_grid, temp_grid.Ny, temp_grid.Nx, temp_grid.Nz, zs[0], ze[0], indx, x_afoil, y_afoil, flat_portion, semicircle_portion, x_circle_start, len_y, delta_y, n_iterations);
    auto rot_x = Create3DMatrix(temp_grid.Nx, temp_grid.Ny, temp_grid.Nz); auto rot_y = Create3DMatrix(temp_grid.Nx, temp_grid.Ny, temp_grid.Nz); auto rot_z = Create3DMatrix(temp_grid.Nx, temp_grid.Ny, temp_grid.Nz);

    {
        const double y_angle = RotateAngle * (3.14159265358979323846 / 180.0); 
        const double z_angle = 0.0;                                            

        const double cos_y = std::cos(y_angle);
        const double sin_y = std::sin(y_angle);
        const double cos_z = std::cos(z_angle);
        const double sin_z = std::sin(z_angle);

        for (int k = 0; k < Nk; ++k) {
            for (int j = 0; j <Nj; ++j) {
                for (int i = 0; i < Ni; ++i) {
                    const double x_orig = temp_grid.x[i][j][k];
                    const double y_orig = temp_grid.y[i][j][k];
                    const double z_orig = temp_grid.z[i][j][k];

                    //about Z
                    const double x1 =  x_orig * cos_z - y_orig * sin_z;
                    const double y1 =  x_orig * sin_z + y_orig * cos_z;
                    const double z1 =  z_orig;

                    //Y
                    const double x2 =  x1 * cos_y + z1 * sin_y;
                    const double y2 =  y1;
                    const double z2 = -x1 * sin_y + z1 * cos_y;

                    rot_x[i][j][k] = x2;
                    rot_y[i][j][k] = y2;
                    rot_z[i][j][k] = z2;
                }
            }
        }
    }
/////end of grid generation for wing grid/////
    //Block 1 =cap grid
    std::ofstream fout("finite_wing.xyz");
    if (!fout) {
        std::cerr << "Error: cannot open output file: finite_wing.xyz\n";
        return 1;
    }

    //Number of blocks
    fout << "2\n";

    //Dimensions for each block
    fout << Nr << " " << Ntheta << " " << newNz << "\n";
    fout << Ni << " " << Nj << " " << Nk << "\n";

// Block 1 =cap grid
//X
for(int k=0;k<newNz;++k)
for(int j=0;j<Ntheta;++j)
for(int i=0;i<Nr;++i)
fout << setprecision(10) << X[k][j][i] << "\n";

//Y
for(int k=0;k<newNz;++k)
for(int j=0;j<Ntheta;++j)
for(int i=0;i<Nr;++i)
fout << setprecision(10) << Y[k][j][i] << "\n";

//Z
for(int k=0;k<newNz;++k)
for(int j=0;j<Ntheta;++j)
for(int i=0;i<Nr;++i)
fout << setprecision(10) << Z[k][j][i] << "\n";

//Block 2 =wing grid
    //X
for (int k = 0; k < Nk; ++k) 
for (int j = 0; j <Nj; ++j) 
for (int i = 0; i < Ni; ++i) 
fout << std::setprecision(10) << rot_x[i][j][k] << "\n";
    //Y
for (int k = 0; k < Nk; ++k) 
for (int j = 0; j <Nj; ++j) 
for (int i = 0; i < Ni; ++i)
fout << std::setprecision(10) << rot_y[i][j][k] << "\n";
    //Z
for (int k = 0; k < Nk; ++k) 
for (int j = 0; j <Nj; ++j) 
for (int i = 0; i < Ni; ++i)
fout << std::setprecision(10) << rot_z[i][j][k] << "\n";

    fout.close();
    std::cout << "finite_wing.xyz\n";

    return 0;
}

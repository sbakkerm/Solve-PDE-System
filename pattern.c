/* LIBRARIES */
#include <stdio.h>
#include <math.h>
#include <string.h>

/* CONSTANTS */
#define N 2601 // Number of elements in u,v arrays (n_x^2)
#define n_x 51 // Grid size in x and y dimensions
#define n_t 10000 // Number of time steps
double f = 0.05; // "Feed" rate
double k = 0.063; // "Kill" rate
double r_u = 0.21; // Diffusion coefficient for u
double r_v = 0.105; // Diffusion coefficient for v
double delta_x = 1.6; // Grid spacing
double delta_t = 0.5; // Time spacing

/* FUNCTIONS */
void PrintSolution(double v[]);
void InitializeArrays(double u[], double v[]);
double ComputeDifferentialReflective(double arr[], int i, int j);
double Update_u(double u[], double v[], int i, int j, int isPeriodic);
double Update_v(double u[], double v[], int i, int j, int isPeriodic);
void IterativeSolution(double u[], double v[], int isPeriodic);
double ComputeDifferentialPeriodic(double arr[], int i, int j);

//int argc, char *argv []
/* MAIN */
int main(int argc, char *argv []){

	/*
	 * Check command arguments
	 */
	int isPeriodic = 0; // False

	if (argc > 1){

		char* str = argv[1];

		// If "p" is an input argument
		if (strcmp(str, "p") == 0){

			isPeriodic = 1; // True
		}
	}

	/*
	 * Initialize arrays
	 */
	double u[N];
	double v[N];
	InitializeArrays(u, v);

	// Initialize number of time steps undergone
	int timeStep = 0;

	// Perform n_t time steps
	while (timeStep < n_t){

		// Solve for u, v
		IterativeSolution(u, v, isPeriodic);

		// Increment the number of time steps undergone
		timeStep += 1;
	}

	// Print final v
	printf("\n");
	PrintSolution(v);

	return 0;
}

/*
 * Prints a 2D representation of the input 1D array v.
 */
void PrintSolution(double v[]){

	// Iterate through rows of 2D image
	for (int i = 0; i < n_x; i++){

		// Iterate through columns of 2D image
		for (int j = 0; j < n_x; j++){

			// If the element's value is greater than 0.25
			if (v[n_x*i + j] > 0.25){

				printf("X");
			}
			// If the element's value is less than or equal to 0.25
			else{

				printf(" ");
			}
		}
	printf("\n");
	}
}

/*
 * Creates the arrays u, v using their initial conditions.
 * Takes the 1D double arrays u, v as input.
 */
void InitializeArrays(double u[], double v[]){

	// Initialize indexing variable
	int index;

	/*
	 * Initialize elements of u, v
	 */
	// Iterate through rows of 2D image
	for (int i = 0; i < n_x; i++){

		// Iterate through columns of 2D image
		for (int j = 0; j < n_x; j++){

			// Compute the index value (i, j)
			index = n_x*i + j;

			// If i,j is in center box (initial condition)
			if ( ((i >= 10) && (i <= 20)) && ((j >= 10) && (j <= 20)) ){

				// Set values of u, v at (i, j)
				u[index] = 0.0;
				v[index] = 1.0;
			}
			// If i,j is anywhere else
			else{

				// Set values of u, v at (i, j)
				u[index] = 1.0;
				v[index] = 0.0;
			}
		}
	}

}



/*
 * Computes the differential (d^2 u/dx^2 + d^2 u/dy^2) OR (d^2 v/dx^2 + d^2 v/dy^2).
 * Takes a double 1D array as input (can be 'u' or 'v') and the 2D indices i and j as input.
 * Uses Reflective Boundary conditions.
 */
double ComputeDifferentialReflective(double arr[], int i, int j){

	// Initialize the quantities for computing the differential
	double east, west, north, south;
	double center = arr[n_x*i + j];

	/*
	 * If i and/or j is at the boundary, then the quantities are adjusted
	 * to agree with a reflective boundary.
	 */
	if (i == 0){
		east = 2*arr[n_x*(i+1) + j];
		west = 0.0;
	}
	else if (i == (n_x - 1)){
		east = 0.0;
		west = 2*arr[n_x*(i-1) + j];
	}
	else {
		east = arr[n_x*(i+1) + j];
		west = arr[n_x*(i-1) + j];
	}

	if (j == 0){
		north = 2*arr[n_x*i + (j+1)];
		south = 0.0;
	}
	else if (j == (n_x - 1)){
		north = 0.0;
		south = 2*arr[n_x*i + (j-1)];
	}
	else{
		north = arr[n_x*i + (j+1)];
		south = arr[n_x*i + (j-1)];
	}

	// Compute the differential (d^2 u/dx^2 + d^2 u/dy^2) OR (d^2 v/dx^2 + d^2 v/dy^2)
	double differential = (east + west + north + south - 4*center)/pow(delta_x, 2);

	return differential;
}

/*
 * Computes the differential (d^2 u/dx^2 + d^2 u/dy^2) OR (d^2 v/dx^2 + d^2 v/dy^2).
 * Takes a double 1D array as input (can be 'u' or 'v') and the 2D indices i and j as input.
 * Uses Periodic Boundary conditions.
 */
double ComputeDifferentialPeriodic(double arr[], int i, int j){

	// Initialize the quantities for computing the differential
	double east, west, north, south;
	double center = arr[n_x*i + j];

	/*
	 * If i and/or j is at the boundary, then the quantities are adjusted
	 * to agree with a periodic boundary.
	 */
	if (i == 0){
		east = arr[n_x*(i+1) + j];
		west = arr[n_x*(n_x - 1) + j];
	}
	else if (i == (n_x - 1)){
		east = arr[j];
		west = arr[n_x*(i-1) + j];
	}
	else {
		east = arr[n_x*(i+1) + j];
		west = arr[n_x*(i-1) + j];
	}

	if (j == 0){
		north = arr[n_x*i + (j+1)];
		south = arr[n_x*i + (n_x - 1)];
	}
	else if (j == (n_x - 1)){
		north = arr[n_x*i];
		south = arr[n_x*i + (j-1)];
	}
	else{
		north = arr[n_x*i + (j+1)];
		south = arr[n_x*i + (j-1)];
	}

	// Compute the differential (d^2 u/dx^2 + d^2 u/dy^2) OR (d^2 v/dx^2 + d^2 v/dy^2)
	double differential = (east + west + north + south - 4*center)/pow(delta_x, 2);

	return differential;
}





/*
 * Computes and returns the updated double value of 'u' at the 2D index (i, j).
 * Takes 1D double arrays u and v as input and 2D indices i and j as input.
 */
double Update_u(double u[], double v[], int i, int j, int isPeriodic){

	// Compute the differential for updating u(i,j)
	double differential;
	if (isPeriodic == 1){
		differential = ComputeDifferentialPeriodic(u, i, j);
	}
	else {
		differential = ComputeDifferentialReflective(u, i, j);
	}

	// Compute the 1D index from the 2D indices i and j
	int index = n_x*i + j;

	// Compute the updated double value u(i,j) represented by the variable 'quantity'
	double quantity = f*(1 - u[index]) - u[index]*pow(v[index], 2) + r_u*differential;
	quantity *= delta_t;
	quantity += u[index];

	return quantity;
}

/*
 * Computes and returns the updated  double value of 'v' at the 2D index (i, j).
 * Takes 1D double arrays u and v as input and 2D indices i and j as input.
 */
double Update_v(double u[], double v[], int i, int j, int isPeriodic){

	// Compute the differential for updating v(i,j)
	double differential;
	if (isPeriodic == 1){
		differential = ComputeDifferentialPeriodic(v, i, j);
	}
	else {
		differential = ComputeDifferentialReflective(v, i, j);
	}

	// Compute the 1D index from the 2D indices i and j
	int index = n_x*i + j;

	// Compute the updated double value v(i,j) represented by the variable 'quantity'
	double quantity = u[index]*pow(v[index], 2) - (f + k)*v[index] + r_v*differential;
	quantity *= delta_t;
	quantity += v[index];

	return quantity;
}

/*
 * Iterates through all the elements in u and v and updates their values.
 * Takes 1D double arrays u and v as input.
 */
void IterativeSolution(double u[], double v[], int isPeriodic){

	// Initialize the variables to store the updated values
	double u_updated;
	double v_updated;
	double new_u[n_x*n_x];
	double new_v[n_x*n_x];

	// Iterate through rows of the 2D image
	for (int i = 0; i < n_x; i++){

		// Iterate through columns of the 2D image
		for (int j = 0; j < n_x; j++){

			// Compute and store the updated values at 2D index (i,j)
			u_updated = Update_u(u, v, i, j, isPeriodic);
			v_updated = Update_v(u, v, i, j, isPeriodic);

			// save the updated values for (i,j) to a temp matrix
			new_u[n_x*i + j] = u_updated;
			new_v[n_x*i + j] = v_updated;
		}
	}

	// Iterate through rows of the 2D image
	for (int i = 0; i < n_x; i++){

		// Iterate through columns of the 2D image
		for (int j = 0; j < n_x; j++){

			// Replace the old values at 2D index (i,j) with the updated values
			u[n_x*i + j] = new_u[n_x*i + j];
			v[n_x*i + j] = new_v[n_x*i + j];
		}
	}
}






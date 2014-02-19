
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include "psflib.h"
#include "pdblib.h"
#include "stringlib.h"

using namespace std;

int parse_command_line (int, char **, FILE **, FILE **, FILE **, double *, double *);
void compute_grid_dim(double **, int, double, double, int *, double *);
void compute_esp_grid(double **, int, double *, double , double , double ***, int *, double *);
void print_dx_file(FILE *, double ***, int *, double *, double);

int main (int argc, char* argv[]) {

	int i, j;
	FILE *psfFile; //pointer to psf file
	FILE *pdbFile; //pointer to pdb file
	FILE *dxFile; // pointer to dx file
	double cut=12.0; //cutoff for grid assigned default value of 12.0
	double delta=1.0; //grid spacing assigned a default value of 1.0
	double *charges;
	double **coord;
	double ***grid;
	double origin[3];
	int gridDim[3];
	int nAtoms;
	char buffer[1024];
	FILE *temp1;

	// read command line arguments
	parse_command_line(argc, argv, &psfFile, &pdbFile, &dxFile, &cut, &delta);

	// read psf file header to get number of atoms
	printf("Reading psf file\n");
	read_psf_header(psfFile,&nAtoms);
	printf("Number of atoms in psf file: %6d\n",nAtoms);

	// allocate arrays 
	charges = new double [nAtoms];
	coord = new double* [nAtoms];
	for (i=0;i<nAtoms;i++) {
		coord[i] = new double [3];
	}

	// read charges from psf file
	read_psf_charges(psfFile,nAtoms,charges);
	fclose(psfFile);

	// read pdb
	read_pdb(pdbFile,coord,nAtoms);
	fclose(pdbFile);

	// compute ESP on grid - grid depends on coordinates, delta and cutoff
	compute_grid_dim(coord,nAtoms,cut,delta,gridDim,origin);
	printf("Size of grid: %6d, %6d, %6d\n",gridDim[0],gridDim[1],gridDim[2]);
	printf("Grid origin: %6.3f%6.3f%6.3f\n",origin[0],origin[1],origin[2]);
	grid = new double** [gridDim[0]];
	for (i=0;i<gridDim[0];i++) {
		grid[i] = new double* [gridDim[1]];
		for(j=0;j<gridDim[1];j++) {
			grid[i][j] = new double [gridDim[2]];
		}
	}
	compute_esp_grid(coord,nAtoms,charges,cut,delta,grid,gridDim,origin);
	
	// print dx file
	print_dx_file(dxFile,grid,gridDim,origin,delta);
	fclose(dxFile);	

}

void print_dx_file(FILE *dxFile, double ***grid, int *gridDim, double *origin, double delta) {

	int x, y, z;
	int counter;

	fprintf(dxFile,"# electrostatic potential from psf and pdb file\n");
	fprintf(dxFile,"# computed with program written by Martin McCullagh 2/19/14\n");
	fprintf(dxFile,"# potential is output in kT/e\n");
	fprintf(dxFile,"object 1 class gridpositions counts %6d %6d %6d\n",gridDim[0],gridDim[1],gridDim[2]);
	fprintf(dxFile,"origin %14.6e%14.6e%14.6e\n",origin[0],origin[1],origin[2]);
	fprintf(dxFile,"delta %14.6e%14.6e%14.6e\n",delta,0.0,0.0);
	fprintf(dxFile,"delta %14.6e%14.6e%14.6e\n",0.0,delta,0.0);
	fprintf(dxFile,"delta %14.6e%14.6e%14.6e\n",0.0,0.0,delta);
	fprintf(dxFile,"object 2 class gridconnections counts %6d %6d %6d\n",gridDim[0],gridDim[1],gridDim[2]);
	fprintf(dxFile,"object 3 class array type double rank 0 items%8d data follows\n",gridDim[0]*gridDim[1]*gridDim[2]);

	counter=0;
	for (x=0;x<gridDim[0];x++) {
		for(y=0;y<gridDim[1];y++) {
			for (z=0;z<gridDim[2];z++) {

				fprintf(dxFile,"%14.6e",grid[x][y][z]);
				counter++;
				if (counter==3) {
					fprintf(dxFile,"\n");
					counter = 0;
				}

			}
		}
	}
	// might need new line
	if (counter!=0) {
		fprintf(dxFile,"\n");
	}

	// print end of dx file
	fprintf(dxFile,"attribute \"dep\" string \"positions\"\n");
	fprintf(dxFile,"object \"regular positions regular connections\" class field\n");
	fprintf(dxFile,"component \"positions\" value 1\n");
	fprintf(dxFile,"component \"connections\" value 2\n");
	fprintf(dxFile,"component \"data\" value 3\n");

}

// compute the esp on a grid
void compute_esp_grid(double **coord, int nAtoms, double *charges, double cut, double delta, double ***grid, int *gridDim, double *origin) {

	double ke = 196.748; //kT*Angstroms/e^2 
	double innerCut = 2.0;
	int atom;
	int i, j, k;
	int x, y, z;
	int cutGrid;
	int atomGrid[3];
	int currentGrid[3];
	double dist;
	double temp;

	cutGrid = int(cut/delta);

	//zero grid
	for(x=0;x<gridDim[0];x++) {
		for(y=0;y<gridDim[0];y++) {
			for (z=0;z<gridDim[0];z++) {
				grid[x][y][z]=0;
			}
		}
	}

	printf("Starting atom loop...\n");

	// loop through all atoms and compute their contributions to surrounding grid points
	for (atom=0;atom<nAtoms;atom++) {
		for (j=0;j<3;j++) {
			atomGrid[j] = int((coord[atom][j]-origin[j])/delta);
		}
		
		// Move forward along grid
		for (x=-cutGrid;x<cutGrid;x++) {
			currentGrid[0] = atomGrid[0] + x;
			if (currentGrid[0]>=0 && currentGrid[0]<gridDim[0]) {
				for (y=-cutGrid;y<cutGrid;y++) {
					currentGrid[1] = atomGrid[1] + y;
					if (currentGrid[1]>=0 && currentGrid[1]<gridDim[1]) {
						for (z=-cutGrid;z<cutGrid;z++) {
							currentGrid[2] = atomGrid[2] + z;
							if (currentGrid[2]>=0 && currentGrid[2]<gridDim[2]) {
								dist = 0;
								for (j=0;j<3;j++) {
									temp = currentGrid[j]*delta+origin[j]-coord[atom][j];
									dist += temp*temp;
								}
								dist = sqrt(dist);
								if (dist>innerCut) {
									grid[currentGrid[0]][currentGrid[1]][currentGrid[2]] += ke*charges[atom]/dist;
								} else {
//									printf("grid[%3d][%3d][%3d] = %10.5f\n",currentGrid[0],currentGrid[1],currentGrid[2],ke*charges[atom]/dist);
								}
							}
						}
					}
				}
			}
		} // end x grid loop
	
	} // end atom loop

} // end subroutine compute_esp_grid


// compute the dimensions of the grid based on atomic coordinates, 
void compute_grid_dim(double **coord, int nAtoms, double cut, double delta, int *gridDim, double *origin) {

	int atom;
	int i, j;
	double max[3];
	double min[3];

	for(atom=0;atom<nAtoms;atom++) {
		for (j=0;j<3;j++) {
			if (atom==0 || coord[atom][j]>max[j]) {
				max[j] = coord[atom][j];
			}	
			if (atom==0 || coord[atom][j]<min[j]) {
				min[j] = coord[atom][j];
			}	
		}
	}
	// compute grid dimensions
	for (j=0;j<3;j++) {
		gridDim[j] = int( (max[j]-min[j]+2*cut)/delta)+1;
		origin[j] = min[j]-cut;
	}

}


int parse_command_line (int argc, char* argv[], FILE **psfFile, FILE **pdbFile, FILE **dxFile, double *cut, double *delta) {

	char psfFileName[1024]={'0'};
	char pdbFileName[1024]={'0'};
	char dxFileName[1024]={'0'};
	char buffer[1024];
	int i;

	// First make sure there are a reasonable number of command line options
	if (argc<3) {
		printf("Insufficient command line options\n");
		printf("Usage: psf2dx -psf [input psf file name] -pdb [input pdb file name] -dx [output dx file name] -cut [optional cutoff value] -delta [optional grid spacing value]\n");
		exit(1);
	}

	// Read command line arguments
	for (i=1;i<argc;i++) {
		// if psf file get argument and open file
		if (strcmp(argv[i], "-psf") ==0) {
			strcpy(psfFileName,argv[++i]);
			*psfFile = fopen(psfFileName,"r");
		// if pdb file get argument and open file
		} else if (strcmp(argv[i], "-pdb") ==0) {
			strcpy(pdbFileName, argv[++i]);
			*pdbFile = fopen(pdbFileName,"r");
		// if dx file get argument and open file
		} else if (strcmp(argv[i], "-dx") ==0) {
			strcpy(dxFileName, argv[++i]);
			*dxFile = fopen(dxFileName,"w");
		// if cut read in cutoff value
		} else if (strcmp(argv[i], "-cut") ==0) {
			*cut = atof(argv[++i]);
		// if delta read in delta (grid spacing) value
		} else if (strcmp(argv[i], "-delta") ==0) {
			*delta = atof(argv[++i]);
		} else {
			printf("Unrecognized command line option: %s\n",argv[i]);
			printf("Usage: psf2dx -psf [input psf file name] -pdb [input pdb file name] -dx [output dx file name] -cut [optional cutoff value] -delta [optional grid spacing value]\n");
			exit(1);
		}
	}
	// Check everything is defined and print log file
	if (strncmp(psfFileName,"0",1)==0) {
		printf("psf file is not defined.  Must be defined by -psf [psf File Name]\n");
		exit(1);
	} else {
	       printf("using psf file: %s\n",psfFileName);
	}	       
	if (strncmp(pdbFileName,"0",1)==0) {
		printf("pdb file is not defined.  Must be defined by -pdb [pdb File Name]\n");
		exit(1);
	} else {
	       printf("using pdb file: %s\n",pdbFileName);
	} 
	if (strncmp(dxFileName,"0",1)==0) {
		printf("dx file is not defined.  Must be defined by -dx [dx File Name]\n");
		exit(1);
	} else {
	       printf("using dx file: %s\n",dxFileName);
	}
       printf("Cutoff value: %f\n",*cut);	
       printf("Grid spacing: %f\n",*delta);
		

}


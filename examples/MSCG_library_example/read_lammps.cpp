//read_lammps.cpp
#include "read_lammps.h"

// Internal function prototypes
void rread_lammps_header(rLammpsData* lammps_data, int* const current_n_sites, int *const timestep, double *const time, double* simulation_box_half_lengths, const int dynamic_types, const int dynamic_state_sampling);
int rread_lammps_body(rLammpsData* lammps_data, const int current_n_sites, double* x, double* f, const int dynamic_types, const int dynamic_state_sampling);

// C++ tokenizing function for strings.

int rStringSplit(std::string source, const char *const delimiter, std::string* results)
{
	int count = 0;
    size_t prev = 0;
    size_t next = 0;

    while ((next = source.find_first_of(delimiter, prev)) != std::string::npos)
    {
        if (next - prev != 0)
        {
            results[count] = source.substr(prev, next - prev);
	        count++;
        }
        prev = next + 1;
    }

    if (prev < source.size())
    {
        results[count] = source.substr(prev);
        count++;
    }
    return count;
}

void rfinish_lammps_reading(rLammpsData* lammps_data)
{
    //close trajectory file
    lammps_data->trajectory_stream.close();
    
    //cleanup allocated memory
    delete [] lammps_data->elements;
	delete lammps_data;
}

// Read the initial frame of a lammps-dump-format trajectory.

rLammpsData* rread_initial_lammps_frame(rLammpsData* lammps_data, const int n_cg_sites, int* cg_site_types, double* simulation_box_half_lengths, double* x, double* f, char* filename)
{
	assert(n_cg_sites > 0);	
	lammps_data = new rLammpsData;
	lammps_data->header_size = 0;
	int n_sites = 0;
	int timestep;
	double time;
	
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
	lammps_data->trajectory_stream.open(filename, std::ifstream::in);
	if (lammps_data->trajectory_stream.fail()) {
		printf("Problem opening lammps trajcetory %s\n", filename);
		exit(EXIT_FAILURE);
	}
	
	//read header for first frame 
	rread_lammps_header(lammps_data, &n_sites, &timestep, &time, simulation_box_half_lengths, 0, 0);
	
	if(n_sites <= 0) {
		exit(EXIT_FAILURE);
    }
    
    //allocate position and force vectors
    lammps_data->elements = new std::string[lammps_data->header_size];
    lammps_data->state_pos = -1;
    lammps_data->type_pos = -1;
        
    // Check that the trajectory is consistent with the desired CG model.
    printf("Number of CG sites defined in top.in (%d) is consistent with trajectory.\n", n_cg_sites);
     
    //read the body of the frame into memory
    if ( rread_lammps_body(lammps_data, n_cg_sites, x, f, 0, 0) != 1 ) {
    	printf("Cannot read the first frame!\n");			
    	delete [] lammps_data->elements;			
		exit(EXIT_FAILURE);
    }
 
    return lammps_data;
}

rLammpsData* rread_next_lammps_frame(rLammpsData* lammps_data, int const reference_atoms, double* simulation_box_half_lengths, double* x, double* f)
{
	int return_value = 1;
	int current_n_sites = reference_atoms;
	int current_timestep;
	double time;

	rread_lammps_header(lammps_data, &current_n_sites, &current_timestep, &time, simulation_box_half_lengths, 0, 0);    

 	if (reference_atoms != current_n_sites) {
 		printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
 		return_value = 0;
 	} else if ( rread_lammps_body(lammps_data, current_n_sites, x, f, 0, 0) != 1) {
    	printf("Cannot read the frame at time %lf!\n", time);
    	return_value = 0;
    }

 	// Return 1 if successful, 0 otherwise.
 	return lammps_data;
}

void rread_lammps_header(rLammpsData * lammps_data, int* const current_n_sites, int *const timestep, double *const time, double* simulation_box_half_lengths, const int dynamic_types, const int dynamic_state_sampling)
{
	double low = 0.0;
	double high = 0.0;
	std::string line;
	int flag = 1; 
	
	while(flag == 1) {
		//read next line of header (and wrap-up if end-of-file)
		std::getline(lammps_data->trajectory_stream, line);
		
		//test if it is a labeled line (all LAMMPS labels start with "ITEM:")
		if( line.compare(0, 5, "ITEM:") == 0 ) {
			//find out which label matched (skip space after ITEM:)
			if( line.compare(6, 15, "NUMBER OF ATOMS") == 0) {
				
				//read number of atoms
				lammps_data->trajectory_stream >> *current_n_sites;
				
			} else if( line.compare(6, 10, "BOX BOUNDS") == 0) {
					
				//read in bounds (low high) for 3 dimensions
				for(int pos=0; pos<3; pos++) {
					lammps_data->trajectory_stream >> low >> high;
					simulation_box_half_lengths[pos] = (high - low)/2.0;
				}	
				
			} else if( line.compare(6, 8, "TIMESTEP") == 0) {
				
				//read in timestep value
				lammps_data->trajectory_stream >> *time;
				(*timestep)++;
			
			} else if( line.compare(6, 5, "ATOMS") == 0) {
				
				//read labels for body of frame
				flag = 0; 
				size_t prev = 11;
				size_t next = 0;
				int set_x = 0;
				int set_f = 0;
				int set_type = 0;
				int set_state = 0;
				lammps_data->header_size = 0;
				
				//check for xpos and fpos as we tokenize string to determine number of columns in body
				while ((next = line.find_first_of(" ", prev)) != std::string::npos) {
					if( (next - prev) == 0 ) { //check if empty
						prev++;
						continue;
					} else if( line.compare(prev, 1, "x") == 0 ) {
						lammps_data->x_pos = lammps_data->header_size;
						set_x = 1;
					} else if( line.compare(prev, 2, "fx") == 0 ) {
						lammps_data->f_pos = lammps_data->header_size;
						set_f = 1;
					} else if( line.compare(prev, 4, "type") == 0 ) {		
						lammps_data->type_pos = lammps_data->header_size;
						set_type = 1;
					} else if( line.compare(prev, 5, "state") == 0 ) {
						lammps_data->state_pos = lammps_data->header_size;
						set_state = 1;
					}
					lammps_data->header_size++;
					prev = next;
				}
		
				if (prev < line.size()) {
                	if( line.compare(prev, 1, "x") == 0 ) {
                        lammps_data->x_pos = lammps_data->header_size;
                    	set_x = 1;
                    } else if( line.compare(prev, 2, "fx") == 0 ) {
                        lammps_data->f_pos = lammps_data->header_size;
                    	set_f = 1;
                    } else if( line.compare(prev, 4, "type") == 0 ) {
                        lammps_data->type_pos = lammps_data->header_size;
                    	set_type = 1;
					} else if( line.compare(prev, 5, "state") == 0 ) {
						lammps_data->state_pos = lammps_data->header_size;
						set_state = 1;
					}
                	lammps_data->header_size++;
                }
				
				//verify that necessary information was extracted to input
				if( (set_x == 0) || (set_f == 0) ) {
					printf("Warning: Was not able to find either x position or f position when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
				if ( (dynamic_types == 1) && (set_type == 0) ) {
					printf("Warning: Type information not detected in header when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
				if ( (dynamic_state_sampling == 1) && (set_state == 0) ) {
					printf("Warning: State probability information not detected in header when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
				
			} else {
				printf("Unrecognized line in frame header: %s", line.c_str() );
				
			}	//close inner if/else if structure
		}		//close ITEM match
	}			//close while	
	return;
}

int rread_lammps_body(rLammpsData *const lammps_data, const int current_n_sites, double* x, double* f, const int dynamic_types, const int dynamic_state_sampling)
{
	//read in current_n_sites lines to extract position and force information
	int j = 0;
	int return_value = 1;
	std::string line;
	
	//allocate space for array of strings based on lammps_data->header_size
	for(int i=0; i < current_n_sites; i++)
	{
		//read in next line and tokenize 
		std::getline(lammps_data->trajectory_stream, line, '\n');	
		if(  (j = rStringSplit(line, " \t", lammps_data->elements)) != lammps_data->header_size ) {	//allow for trailing white space
			printf("Warning: Number of fields detected in frame body");
			printf(" (%d) does not agree with number expected from frame header (%d)!\n", j, lammps_data->header_size);
			return_value = -1;
			}
			
		//extract position information
		for(j = 0; j < 3; j++) {
			x[i*3+j] = atof( lammps_data->elements[j + lammps_data->x_pos].c_str() );
		}
		
		//extract force information
		for(j = 0; j < 3; j++) {
			f[i*3+j] = atof( lammps_data->elements[j + lammps_data->f_pos].c_str() );
		}
	}
	return return_value;
}
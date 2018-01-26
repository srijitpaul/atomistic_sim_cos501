#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <stdlib.h>
using namespace std;

double force_i(int i_particle, int bead_index, int dim, double* x);
double etot(double* mom_j);
double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double *mom_out, double *x_out, double *sq_en);
double mass(int index);

double delt         =   0.01;
auto num_particles  =   10;
auto Nd             =   2;
int Nbeads          =   1;
double beta         =   600;
double L            =   8;/*16x16 lattice size */
double near_neigh   =   2.62;/*Angstrom*/
double r_cutoff     =   5.148;/*Angstrom*/
double D0           =   0.316;
double alpha0       =   1.43;/*Angstrom inverse*/
double r0           =   2.34;/*Angstrom*/
/*
//Normal distribution//
double beta = 4;
double sigma = sqrt(2*mass/beta);
std::normal_distribution<double> nd(0,sigma);
*/

/*Calculate mass on ith particle*/

double mass(int index){
    return 1;
}

/*Defining Euclidean norm*/
/*double norm(double* x_init, double* x_final) {

    double norm = 0;
    for(int i = 0; i<Nd; i++){
        norm += pow((x_final[i] - x_init[i]),2);
    }
}
*/
//Generated lattice
void lat_pos(int i, int j, double* R) {

        R[0] = i * near_neigh;
        R[1] = j * near_neigh;
    }

double distance(double* x, double* R){

    double dist = 0;
    for(int i_particle=0; i_particle<num_particles; i_particle++) {
        for(int j=0; j<Nbeads; j++) {
            for (int dim = 0; dim < Nd; dim++) {

                dist += pow(*(x + i_particle * Nbeads * Nd + j * Nd + dim)  - *(R + dim),2);

            }
        }
    }
    dist = sqrt(dist);
    return dist;
}
//Force on the ith particle due to all other particles in one particular dimension
double force_i( int i_particle, int bead_index, int dim, double* x){


    double force = 0;
//    int prev_bead = ((Nbeads + ((bead_index-1)%Nbeads)) % Nbeads);
//    int next_bead = ((Nbeads + ((bead_index+1)%Nbeads)) % Nbeads);

    for(int j = 0; j<num_particles; j++) {

       /* if(j != i_particle)
            force += - mass(j) * (x_j[j]);
        else{}*/

    }
    double R[2];
    double dist;
    double func;
    for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {
                for ( int xcoord = -L; xcoord < L; xcoord ++) {
                
                    for( int ycoord = -L; ycoord < L; ycoord ++){
                        
                        lat_pos(xcoord, ycoord, R);
                        dist = distance(x, R);
                        if(dist >= r_cutoff)
                            dist = 0;
                        else{
                            func = exp(-1 * alpha0 *(dist - r0));
                            //cout<<( func - 1)<<endl; 
                            force += 2 * alpha0 * D0 * func * ( func - 1); 
                        }
                    }
                }
            }
        }
    }
/*
    if(num_particles==1) {

        force = -mass(0) * (pow((1/beta),2)*(2 *(*(x + i_particle * Nbeads * Nd + bead_index * Nd + dim))
                            - (*(x + i_particle * Nbeads * Nd + prev_bead * Nd + dim))
                            - (*(x + i_particle * Nbeads * Nd + next_bead * Nd + dim)))
                            + *(x + i_particle * Nbeads * Nd  + bead_index * Nd + dim)
                            + 0.3 * pow(*(x + i_particle * Nbeads * Nd  + bead_index * Nd + dim), 2)
                            + 0.04 * pow(*(x + i_particle * Nbeads * Nd  + bead_index * Nd + dim), 3));
    }
    */    
    //cout<<"force "<<force<<endl;
    return  force;

}

double etot(double* mom_j){

    double K = 0;
    for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {

                K += pow(*(mom_j + i * Nbeads * Nd + j * Nd + k),2)/(2*mass(i));

            }
        }
    }

    return K;
}
double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double* mom_out, double* x_out, double *sq_en) {

    double avg_en=0;
    double x_dum[num_particles][Nbeads];
    double mom_dum[num_particles][Nbeads];
    double *x   = new double[num_particles * Nbeads * Nd];
    double *mom = new double[num_particles * Nbeads * Nd];
    double *x_temp   = new double[num_particles * Nbeads * Nd];
    double *centroid = new double[num_particles * Nd];



    /*Initial conditions loop*/
    for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
	        x_dum[i][j] = 0;
            mom_dum[i][j] = 0;
            for (int k = 0; k < Nd; k++) {
                *(x + i * Nbeads * Nd  + j * Nd + k) = *(x_init + i * Nbeads * Nd + j * Nd + k);
                *(mom + i * Nbeads * Nd  + j * Nd + k) = *(mom_init + i * Nbeads * Nd + j * Nd + k);
				//cout<<"Momentum bead index "<< j<<"= \t\t"<< *(mom + i * Nbeads * Nd + j * Nd + k)<<endl;

            }
        }
    }
    avg_en += etot(mom);

//    if(constraint==1)constraint_momentum_to_dividing_surface(mom, dxi); 
    for(int i=0; i<=(npart   -   1); i++) {
    	for (int p_index = 0; p_index < num_particles; p_index++) {
        	for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
            	for (int dim = 0; dim < Nd; dim++) {

                	x_dum[p_index][bead_index] = *(x + p_index * Nbeads * Nd + bead_index * Nd + dim) +
                                            (*(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) / mass(p_index)) *
                                            delt +
                                            0.5 * (force_i(p_index, bead_index, dim, x) / mass(p_index))
                                            * pow(delt, 2);
            	}
           	}
		}
    	for (int p_index = 0; p_index < num_particles; p_index++) {
		
			for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
				for (int dim = 0; dim < Nd; dim++) {

					*(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim) = fmod((x_dum[p_index][bead_index] + 10*L),L);
					//cout<<"Momentum from previous step bead index "<<bead_index<<"=\t\t"<< *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;
					//cout<<"Position after verlet step bead index "<<bead_index<<"=\t\t"<< *(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;

				}
			}
		}
		//cout<<"dxi "<<dxi[0]<<endl;
/*		if(constraint == 1 ) {	
			constraint_to_dividing_surface(x_temp, mom, dxi, xi_current);
            constraint_momentum_to_dividing_surface(mom, dxi); 
			get_centroid(x_temp,centroid);
			get_reaction_coordinates(centroid, xi_current, xi, dxi);
		}
*/
    	for (int p_index = 0; p_index < num_particles; p_index++) {
			for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
				for (int dim = 0; dim < Nd; dim++) {

					mom_dum[p_index][bead_index] = *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) +
									  0.5 * (force_i(p_index, bead_index, dim, x) +
											 force_i(p_index, bead_index, dim, x_temp)) / (mass(p_index)) * delt;
					//cout<<"Momentum before constraint bead index "<< bead_index<<"= \t\t"<< mom_dum[p_index][bead_index]<<endl;
					//cout<<"Position after step bead index "<<bead_index<<"- \t\t"<< *(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;
				}
			}
		}
	 	//constraint_momentum_to_dividing_surface(mom, dxi);
 		
    	for (int p_index = 0; p_index < num_particles; p_index++) {
			for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
				for (int dim = 0; dim < Nd; dim++) {


					*(x + p_index * Nbeads * Nd + bead_index * Nd + dim) = *(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim);
					*(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) = mom_dum[p_index][bead_index];
				}
			}
		}
// 		if(constraint == 1) constraint_momentum_to_dividing_surface(mom, dxi); 

        avg_en += etot(mom);
			//cout<<"Energy "<<avg_en<<endl;
        sq_en[0] += pow(etot(mom),2);
      	
		//cout<<"One complete verlet step"<<endl;

	}

	//cout<<"One npart step in verlet"<<endl;
	for (int p_index = 0; p_index < num_particles; p_index++) {
		for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
			for (int dim = 0; dim < Nd; dim++) {
				*(x_out + p_index * Nbeads * Nd + bead_index * Nd + dim) = *(x + p_index * Nbeads * Nd +
																	 bead_index * Nd + dim);
				*(mom_out + p_index * Nbeads * Nd + bead_index * Nd + dim) = *(mom + p_index * Nbeads * Nd +
																	   bead_index * Nd + dim);
				//cout<<"Momentum after one verlet  bead index "<< bead_index<<"= \t\t"<< *(mom_out + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;
			
			}
		}
	}
	//cout<<*(x_out +1)<<endl;
    delete [] x;
    delete [] mom;
    delete [] x_temp;
	delete [] centroid;
    return avg_en;
}          
/*
double RK2_integrate(int npart) {

    double *x, *mom;
    x = new double[npart];
    mom = new double[npart];

    x[0]    =   x_init;
    mom[0]  =   mom_init;

    double kx1, kx2;
    double kv1, kv2;

    for(int i=0; i<(npart   -   1); i++){

        kx1 = mom[i]/mass;
        kv1 = force(x[i])/mass;
        kx2 = kx1 + delt * kv1;
        kv2 = force(x[i] + delt * kx1)/mass;


        x[i+1]      =   x[i] + delt * (kx1 + kx2)/2;
        mom[i+1]    =   mom[i] + mass * delt * (kv2 + kv1)/2;

    }


    double e =  etot(mom[npart - 1], x[npart - 1]);
    delete [] x;
    delete [] mom;

    return e;
}
double RK4_integrate(int npart) {

    double *x,*mom;
    x = new double[npart];
    mom = new double[npart];

    x[0]    =   x_init;
    mom[0]  =   mom_init;

    double kx1, kx2, kx3, kx4;
    double kv1, kv2, kv3, kv4;
    for(int i=0; i<(npart   -   1); i++){

        kx1 = mom[i]/mass;
        kv1 = force(x[i])/mass;
        kx2 = kx1 + delt * kv1/2;
        kv2 = force(x[i] + delt * kx1/2)/mass;
        kx3 = kx1 + delt * kv2/2;
        kv3 = force(x[i] + delt * kx2/2)/mass;
        kx4 = kx1 + delt * kv3;
        kv4 = force(x[i] + delt * kx3);

        x[i+1]      = x[i] + delt * (kx1 + 2 * kx2 + 2 * kx3 + kx4)/6;
        mom[i+1]    = mom[i] + mass * delt * (kv1 + 2 * kv2 + 2 * kv3 + kv4)/6;

    }

    double e =  etot(mom[npart - 1], x[npart - 1]);
    delete [] x;
    delete [] mom;

    return e;
}
*/
int main() {

    cout<<beta<<endl;

    //std::ofstream vel_verlet("vel_verlet_etot.dat", ios::out);
    std::ofstream vel("qm_vel_verlet_fluc.dat", ios::out);
    int thermostat_freq = 100;
    int total_steps = 100;
    double e_vel_verlet = 0;

    /* Initial conditions */
    double* mom_init    =   new double[num_particles * Nbeads * Nd];
    double* x_init      =   new double[num_particles * Nbeads * Nd];


    /*Initial conditions loop*/
   /* for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {
                *(x_init + i * Nbeads * Nd  + j * Nd + k) = 1;
                *(mom_init + i * Nbeads * Nd  + j * Nd + k) = 0;
            }
        }
    }*/


    double* x_out   =   new double[num_particles * Nbeads * Nd];
    double* mom_out =   new double[num_particles * Nbeads * Nd];
    double sq_en[1];
    double sq_e_vel_verlet =0;
    vel<<scientific;
    vel.precision(16);
    int ch=0;
    int num_trajectories = 1;

    double x[5000][2] = {0};

    double corr[5000] = {0};
//Normal distribution//
    for(int traj = 1; traj<=num_trajectories;traj++) {

        /*Initial conditions loop*/
        for(int i=0; i<num_particles; i++) {
            for(int j=0; j<Nbeads; j++) {
                for (int k = 0; k < Nd; k++) {
                    *(x_init + i * Nbeads * Nd  + j * Nd + k) = 1;
                    *(mom_init + i * Nbeads * Nd  + j * Nd + k) = 0;
                }
            }
        }


        cout<<"Starting new trajectory"<<endl;
        e_vel_verlet = 0;
        sq_e_vel_verlet = 0;
        double sigma[num_particles];

        //Thermalization step
        for (int steps = 1; steps <= total_steps; steps++) {

/*        std::cout << "VELOCITY-VERLET: The energy after nsteps:" << steps - 1 << " is "
                  << velocity_verlet_integrate(steps) << std::endl;
        std::cout << "RUNGE-KUTTA - 2: The energy after nsteps:" << steps - 1 << " is " << RK2_integrate(steps) << std::endl;
        std::cout << "RUNGE-KUTTA - 4: The energy after nsteps:" << steps - 1 << " is " << RK4_integrate(steps) << std::endl;
*/
            for(int i=0; i<num_particles; i++) {

                //sigma[i] = sqrt(mass(i)/beta);
                //std::normal_distribution<double> nd(0, sigma[i]);

                for(int j=0; j<Nbeads; j++) {

                    sigma[i] = sqrt(mass(i)/beta);
                    std::normal_distribution<double> nd(0, sigma[i]);

                    for (int k = 0; k < Nd; k++) {
                        std::default_random_engine de(traj*steps*(j+10)*(k+100));//seed
                        *(mom_init + i * Nbeads * Nd  + j * Nd + k) = nd(de);

                    }
                }
            }

            e_vel_verlet += velocity_verlet_integrate(thermostat_freq, mom_init, x_init, mom_out, x_out, sq_en);
            sq_e_vel_verlet = sq_en[0];
            //vel_verlet << beta<<"\t\t"<<delt * steps * thermostat_freq << "\t\t" << e_vel_verlet / (steps * thermostat_freq)i
            //           << std::endl;

            x_init = x_out;
            mom_init = mom_out;
        }
        cout<<"Energy "<<e_vel_verlet/(total_steps * thermostat_freq)<<endl;
    	cout<<"Time "<<total_steps * thermostat_freq * delt<<endl;


        double dum;
        //calculation step
        //x[0] = 0;

        for(int i=0; i<num_particles; i++) {
            for(int j=0; j<Nbeads; j++) {
                for (int k = 0; k < Nd ; k++) {
                    x[0][k] = pow(*(x_out + i * Nbeads * Nd  + j * Nd + k),2);

                    //if(j==1) cout<<x[0]<<endl;
                }
            }
        }
       
        corr[0] = x[0][0] + x[0][1];

        for(int ch_t=1; ch_t<= 3000; ch_t++) {

            x[ch_t][0] = 0;
            x[ch_t][1] = 0;
	    //cout<<"input "<<*(mom_init + 0)<<endl;
            dum = velocity_verlet_integrate(1, mom_init, x_init, mom_out, x_out, sq_en);
           // cout<<"dum "<<*(x_out + 1)<<endl;
            for(int i=0; i<num_particles; i++) {
                for(int j=0; j<Nbeads; j++) {
                    for (int k = 0; k < Nd; k++) {
                         x[ch_t][k] = pow(*(x_out + i * Nbeads * Nd  + j * Nd + k) - x[0][k],2);
			 //cout<<x[ch_t]<<endl;
                    }
                }
            }
            
            corr[ch_t] += x[ch_t][0] + x[ch_t][1];

            x_init = x_out;
            mom_init = mom_out;

        }
        cout<<"Traj number"<<traj<<"\t\t"<<corr[1]<<"\t\t"<<corr[2]<<endl;

        //vel << pow(1/beta,2)<<"\t\t"<<delt * total_steps * thermostat_freq << "\t\t"
        //          << sq_e_vel_verlet/(total_steps * thermostat_freq) - pow(e_vel_verlet / (total_steps * thermostat_freq),2) <<endl;
    }
    for (int ch_t=1; ch_t <=3000;){

        vel<< ch_t*delt << "\t\t"<< corr[ch_t]/num_trajectories<<endl;
        ch_t += 1;
    }

    //system("gnuplot plot 'vel_verlet_fluc.dat'");

}

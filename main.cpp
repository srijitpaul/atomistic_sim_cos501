#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <stdlib.h>
#include <algorithm> 
using namespace std;

double force_i(int i_particle, int bead_index, int dim, double* x);
double etot(double* mom_j);
double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double *mom_out, double *x_out, double *sq_en, double* dxi, double xi_current);
double mass(int index);

double delt         =   0.01;
auto num_particles   =   1;
auto Nd              =   1;
int Nbeads          =   4;
double nbeads 		=   4;
double beta          =   0.25;
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
//Force on the ith particle due to all other particles in one particular dimension
double force_i( int i_particle, int bead_index, int dim, double* x){

    double force = 0;
    int prev_bead = ((Nbeads + ((bead_index-1)%Nbeads)) % Nbeads);
    int next_bead = ((Nbeads + ((bead_index+1)%Nbeads)) % Nbeads);

   /* for(int j = 0; j<num_particles; j++) {

        if(j != i_particle)
            force += - mass(j) * (x_j[j]);
        else{}
    }*/
    if(num_particles==1) {

        force = -mass(0) * (pow((1/beta),2)*(2 *(*(x + i_particle * Nbeads * Nd + bead_index * Nd + dim))
                            - (*(x + i_particle * Nbeads * Nd + prev_bead * Nd + dim))
                            - (*(x + i_particle * Nbeads * Nd + next_bead * Nd + dim)))
                            + *(x + i_particle * Nbeads * Nd  + bead_index * Nd + dim)
                            + 0.3 * pow(*(x + i_particle * Nbeads * Nd  + bead_index * Nd + dim), 2)
                            + 0.04 * pow(*(x + i_particle * Nbeads * Nd  + bead_index * Nd + dim), 3));
    }
    //cout<<force<<endl;
    return  force;

}

void get_reaction_coordinates(double* qctemp, double xi_current, double* xi_new, double* dxi_new){

	for(int i=0; i<num_particles; i++) {
		for(int k=0; k<Nd; k++) {
			
			xi_new[0]  =  *(qctemp + i * Nd + k) - xi_current;
			*(dxi_new + i * Nd + k) =  1/nbeads;
		}
	}
	//cout<<"xi new "<<xi_new[0]<<endl;
	//if(xi_new[0] >1e8) std::terminate();
		
}
void get_centroid(double* x, double* centroid) {

    for(int i=0; i<num_particles; i++) {
        for(int k=0; k<Nd; k++) {
            for (int j = 0; j < Nbeads; j++) {
				*(centroid + i * Nd + k) += *(x + i * Nbeads * Nd + j * Nd + k)/nbeads;
				//cout<<"centroid particle num "<<i<<" dim "<<k<<" num beads "<<j<<" = \t\t"<<*(centroid + i * Nd + k)<<endl;
	    	}
		}
    }
} 

void constraint_momentum_to_dividing_surface(double* mom, double* dxi) {

	double coeff1 = 0, coeff2 = 0, lambda;
    for(int i=0; i<num_particles; i++) {
        for(int k=0; k<Nd; k++) {

			coeff2 += *(dxi + i * Nd + k) * (*(dxi + i * Nd + k))/mass(i);
			for (int j = 0; j < Nbeads; j++) {
				
				coeff1 += *(dxi + i * Nd + k) * (*(mom + i * Nbeads * Nd + j * Nd + k))/mass(i);
			}
		}
    }
	
	lambda = -coeff1/coeff2/nbeads;
    
	for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {

            	*(mom + i * Nbeads * Nd + j * Nd + k) += lambda * (*(dxi + i * Nd + k));
				//cout<<"Momentum after constraint bead index "<< j<<"= \t\t"<< *(mom + i * Nbeads * Nd + j * Nd + k)<<endl;

            }
        }
    }
}


void constraint_to_dividing_surface(double* x, double* mom, double* dxi, double xi_current) {

    double *centroid = new double[num_particles * Nd];
    std::fill_n(centroid, num_particles * Nd, 0.0f);
    double *qctemp   = new double[num_particles * Nd];
    std::fill_n(qctemp, num_particles * Nd, 0);
    double xi_new[1], dx, dsigma;
	xi_new[0] = 0;
    double *dxi_new  = new double[num_particles * Nd];
    double *d2xi_new = new double[num_particles * Nd * num_particles * Nd];
    double coeff, sigma; 
	//cout<<"called 1"<<endl;
    get_centroid(x,centroid);
    
    //Langrange Multiplier
    double mult = 0;

    int maxiter = 100000;
	int ch = 0;
    for(int iter = 1; iter<=maxiter && ch==0; iter++ ) {

		coeff = mult * delt * delt/nbeads;
		for(int i=0; i<num_particles; i++) {
			for(int k=0; k<Nd; k++) {
						
				*(qctemp + i * Nd + k) = *(centroid + i * Nd + k) + coeff * (*(dxi + i * Nd + k))/mass(i);
				//cout<<"test "<< coeff*(*(dxi + i * Nd + k))/mass(i)<<endl;
				//if(coeff!= coeff)std::terminate();
			}
		}
					 
		get_reaction_coordinates(qctemp, xi_current, xi_new, dxi_new);
		sigma = xi_new[0];
		dsigma = 0;

		for(int i=0; i<num_particles; i++) {
			for(int k=0; k<Nd; k++) {
						
				dsigma += *(dxi_new + i * Nd + k) * delt * delt * (*(dxi + i * Nd + k))/(mass(i)* nbeads);
				//cout<<"dxi_new "<<*(dxi_new + i * Nd + k)<<endl;
				
			}
		}
		
		dx = sigma/dsigma;
		//cout<<"sigma "<<sigma<<endl;
        mult -= dx;
		
		if(fabs(dx)<pow(1,-8) || fabs(sigma) < pow(1,-10)) ch=1;

		if(iter == maxiter) {
			cout<<"Error: SHAKE exceeded maximum number of iterations. Restarting trajectory."<<endl;
			cout<<"dx = "<<dsigma<<"\t\t"<<"sigma = "<<sigma<<endl;
		}
	
	}
    for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {

            	*(x + i * Nbeads * Nd + j * Nd + k) += coeff/mass(i) * (*(dxi + i * Nd + k));
            	*(mom + i * Nbeads * Nd + j * Nd + k) += mult *(delt/nbeads) * (*(dxi + i * Nd + k));
				// cout<<"dxi bead index "<< j<<"= \t\t"<<  mult * (delt/Nbeads) <<endl;

            }
        }
    }
	
	delete [] qctemp;
	delete [] centroid;
	delete [] dxi_new;
	delete [] d2xi_new;
	 
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

double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double* mom_out, double* x_out, double *sq_en, double* dxi, double* xi, double xi_current) {

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

    constraint_momentum_to_dividing_surface(mom, dxi); 
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

					*(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim) = x_dum[p_index][bead_index];
					//cout<<"Momentum from previous step bead index "<<bead_index<<"=\t\t"<< *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;
					//cout<<"Position after verlet step bead index "<<bead_index<<"=\t\t"<< *(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;

				}
			}
		}
		//cout<<"dxi "<<dxi[0]<<endl;	
		constraint_to_dividing_surface(x_temp, mom, dxi, xi_current);
		get_centroid(x_temp,centroid);
		get_reaction_coordinates(centroid, xi_current, xi, dxi);


    	for (int p_index = 0; p_index < num_particles; p_index++) {
			for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
				for (int dim = 0; dim < Nd; dim++) {

					mom_dum[p_index][bead_index] = *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) +
									  0.5 * (force_i(p_index, bead_index, dim, x) +
											 force_i(p_index, bead_index, dim, x_temp)) / (mass(p_index)) * delt;
					//cout<<"Momentum before constraint bead index "<< bead_index<<"= \t\t"<< mom_dum[bead_index]<<endl;
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
 		constraint_momentum_to_dividing_surface(mom, dxi); 

    	for (int p_index = 0; p_index < num_particles; p_index++) {
			for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
				for (int dim = 0; dim < Nd; dim++) {

					//cout<<"Momentum after constraint bead index "<< bead_index<<"= \t\t"<< *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;
					//cout<<"Position after constraint bead index "<<bead_index<<"= \t\t"<< *(x_temp + p_index * Nbeads * Nd + bead_index * Nd + dim)<<endl;
				}
			}
		}
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

    std::ofstream vel_verlet("vel_verlet_etot.dat", ios::out);
    std::ofstream vel("constrainted_qm_vel_verlet.dat", ios::out);
    int thermostat_freq = 100;
    int total_steps = 100;
    double e_vel_verlet = 0;
    double *dxi  = new double[num_particles * Nd];
    std::fill_n(dxi, num_particles * Nd, 0.25f);
	cout<<"dxi "<<dxi[0]<<endl;
	double xi_current = 0.5;
	double xi[1];
    /* Initial conditions */
    double* mom_init    = new double[num_particles * Nbeads * Nd];
    double* x_init      = new double[num_particles * Nbeads * Nd];


    /*Initial conditions loop*/
   /* for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {
                *(x_init + i * Nbeads * Nd  + j * Nd + k) = 1;
                *(mom_init + i * Nbeads * Nd  + j * Nd + k) = 0;
            }
        }
    }*/


    double* x_out   = new double[num_particles * Nbeads * Nd];
    double* mom_out = new double[num_particles * Nbeads * Nd];
    double sq_en[1];
    double sq_e_vel_verlet =0;
    vel<<scientific;
    vel.precision(16);
    int ch=0;
    int num_trajectories = 2;

    double x[3000] = {0};

    double corr[3000] = {0};
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

//        std::cout << "VELOCITY-VERLET: The energy after nsteps:" << steps - 1 << " is "
//                  << velocity_verlet_integrate(steps) << std::endl;
//        std::cout << "RUNGE-KUTTA - 2: The energy after nsteps:" << steps - 1 << " is " << RK2_integrate(steps) << std::endl;
//        std::cout << "RUNGE-KUTTA - 4: The energy after nsteps:" << steps - 1 << " is " << RK4_integrate(steps) << std::endl;
//
            for(int i=0; i<num_particles; i++) {

                //sigma[i] = sqrt(mass(i)/beta);
                //std::normal_distribution<double> nd(0, sigma[i]);

                for(int j=0; j<Nbeads; j++) {

                    sigma[i] = sqrt(mass(i)/beta);
                    std::normal_distribution<double> nd(0, sigma[i]);

                    for (int k = 0; k < Nd; k++) {
                        std::default_random_engine de(traj*steps*(j+10));//seed
                        *(mom_init + i * Nbeads * Nd  + j * Nd + k) = nd(de);

                    }
                }
            }
			//cout<<"xi "<<xi_current<<endl;
            e_vel_verlet += velocity_verlet_integrate(thermostat_freq, mom_init, x_init, mom_out, x_out, sq_en, dxi, xi,  xi_current);
            sq_e_vel_verlet = sq_en[0];
            //vel_verlet << beta<<"\t\t"<<delt * steps * thermostat_freq << "\t\t" << e_vel_verlet / (steps * thermostat_freq)i
            //           << std::endl;

            //x_init = x_out;
            //mom_init = mom_out;
			//xi_current = xi[0];
            for(int i=0; i<num_particles; i++) {
                for(int j=0; j<Nbeads; j++) {
                    for (int k = 0; k < Nd; k++) {
                        *(x_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_out + i * Nbeads * Nd  + j * Nd + k);
                        *(mom_init + i * Nbeads * Nd  + j * Nd + k) =  *(mom_out + i * Nbeads * Nd  + j * Nd + k);
                    }
                }
            }

        }
        cout<<"Energy "<<e_vel_verlet/(total_steps * thermostat_freq)<<endl;
		cout<<"Time "<<total_steps * thermostat_freq * delt<<endl;


        double dum;
        //calculation step
        x[0] = 0;

        for(int i=0; i<num_particles; i++) {
            for(int j=0; j<Nbeads; j++) {
                for (int k = 0; k < Nd; k++) {
                    x[0] += *(x_out + i * Nbeads * Nd  + j * Nd + k);

                    //if(j==1) cout<<x[0]<<endl;
                }
            }
        }
        corr[0] += x[0] * x[0]/(pow(Nbeads,2));

        for(int ch_t=1; ch_t<= 3000; ch_t++) {

            x[ch_t] = 0;
	    //cout<<"input "<<*(mom_init + 0)<<endl;
            dum = velocity_verlet_integrate(1, mom_init, x_init, mom_out, x_out, sq_en, dxi, xi, xi_current);
           // cout<<"dum "<<*(x_out + 1)<<endl;
            for(int i=0; i<num_particles; i++) {
                for(int j=0; j<Nbeads; j++) {
                    for (int k = 0; k < Nd; k++) {
                         x[ch_t] += *(x_out + i * Nbeads * Nd  + j * Nd + k);
			 //cout<<x[ch_t]<<endl;
                    }
                }
            }
            corr[ch_t] += (x[0] * x[ch_t])/16;
	    	if(corr[ch_t]/traj>10000000000){
	        cout<<x[ch_t]<<endl;
			cout<<*(x_init + 0)<<"\t\t"<<*(x_init + 1)<<"\t\t"<<*(x_init + 2)<<"\t\t"<<*(x_init + 3)<<endl;
			//std::terminate();
	    	}

            //x_init = x_out;
            //mom_init = mom_out;
			//xi_current = xi[0];
            for(int i=0; i<num_particles; i++) {
                for(int j=0; j<Nbeads; j++) {
                    for (int k = 0; k < Nd; k++) {
                        *(x_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_out + i * Nbeads * Nd  + j * Nd + k);
                        *(mom_init + i * Nbeads * Nd  + j * Nd + k) =  *(mom_out + i * Nbeads * Nd  + j * Nd + k);
                    }
                }
            }
            

        }
        cout<<"Traj number"<<traj<<"\t\t"<<corr[0]/traj<<endl;

        //vel << pow(1/beta,2)<<"\t\t"<<delt * total_steps * thermostat_freq << "\t\t"
        //          << sq_e_vel_verlet/(total_steps * thermostat_freq) - pow(e_vel_verlet / (total_steps * thermostat_freq),2) <<endl;
    }
    for (int ch_t=0; ch_t <=3000;){

        vel<< ch_t*delt << "\t\t"<< corr[ch_t]/2000<<endl;
        ch_t += 1;
    }

    //system("gnuplot plot 'vel_verlet_fluc.dat'");
	delete [] dxi;
	delete [] x_init;
	delete [] x_out;
	delete [] mom_init;
	delete [] mom_out;
}

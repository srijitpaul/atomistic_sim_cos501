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
int Nbeads          =   128;
double nbeads 		=   128;
double beta          =   0.03125;
double A            = -18/M_PI;
double B            = 13.5/M_PI;
double a            = 8/sqrt(3*M_PI);
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
                            - 2* B * tanh(*(x + i_particle * Nbeads * Nd + next_bead * Nd + dim)) /(a * cosh(*(x + i_particle * Nbeads * Nd + next_bead * Nd + dim)))
                            + (2 * A * exp ( (-2*(*(x + i_particle * Nbeads * Nd + next_bead * Nd + dim)))/a))/(a * pow(1 +
                                exp ( (-2*(*(x + i_particle * Nbeads * Nd + next_bead * Nd + dim)))/a),2)));
    }
    //cout<<force<<endl;
    return  force;

}

void get_reaction_coordinates(double* qctemp, double xi_current, double* xi_new, double* dxi_new){

	for(int i=0; i<num_particles; i++) {
		for(int k=0; k<Nd; k++) {
			
			xi_new[0]  =  *(qctemp + i * Nd + k) - xi_current;
			*(dxi_new + i * Nd + k) =  1;
		}
	}
	//cout<<"xi new "<<xi_new[0]<<endl;
	//if(xi_new[0] >1e8) std::terminate();
		
}
void get_centroid(double* x, double* centroid) {

    std::fill_n(centroid, num_particles * Nd, 0.0f);
    
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

    int maxiter = 1000000;
	int ch = 0;
    for(int iter = 1; iter<=maxiter && ch==0; iter++ ) {

		coeff = (mult)/nbeads;
		for(int i=0; i<num_particles; i++) {
			for(int k=0; k<Nd; k++) {
						
				*(qctemp + i * Nd + k) = *(centroid + i * Nd + k) + coeff * (*(dxi + i * Nd + k))/mass(i);
				//cout<<"centroid "<<  *(qctemp + i * Nd + k) <<endl;
				if(coeff!= coeff){cout<<coeff<<endl;std::terminate();}
			}
		}
					 
		get_reaction_coordinates(qctemp, xi_current, xi_new, dxi_new);
		sigma = xi_new[0];
		dsigma = 0;

		for(int i=0; i<num_particles; i++) {
			for(int k=0; k<Nd; k++) {
						
				dsigma += *(dxi_new + i * Nd + k) * (*(dxi + i * Nd + k))/(mass(i)* nbeads);
				//cout<<"dsigma "<<dsigma<<endl;
				
			}
		}
		
		dx = sigma/dsigma;
		//cout<<"dx "<<fabs(dx)<<endl;
        //cout<<"sigma "<<fabs(sigma)<<endl;
        mult = mult - dx;
		
		if(abs(dx)<1e-8 || abs(sigma) < 1e-10) {ch=1;}

		if(iter == maxiter) {
			cout<<"Error: SHAKE exceeded maximum number of iterations. Restarting trajectory."<<endl;
			cout<<"dx = "<<dsigma<<"\t\t"<<"sigma = "<<sigma<<endl;
		}
	    //cout<<"Iter number "<<iter<<endl;
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

double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double* mom_out, double* x_out, double *sq_en, double constraint, double* dxi, double* xi, double xi_current) {

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

    if(constraint==1)constraint_momentum_to_dividing_surface(mom, dxi); 
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
		if(constraint == 1 ) {	
			constraint_to_dividing_surface(x_temp, mom, dxi, xi_current);
            constraint_momentum_to_dividing_surface(mom, dxi); 
			get_centroid(x_temp,centroid);
			get_reaction_coordinates(centroid, xi_current, xi, dxi);
		}

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
 		if(constraint == 1) constraint_momentum_to_dividing_surface(mom, dxi); 

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
double heavside(double x) {
	if(x<0) return 0;
	else return 1;
}

int main() {

    cout<<beta<<endl;
	double constraint;
    //std::ofstream vel_verlet("vel_verlet_etot.dat", ios::out);
    std::ofstream vel("neg_1point0_constrainted_qm_vel_verlet.dat", ios::out);
    int thermostat_freq = 100;
    int total_steps = 1000;
    double e_vel_verlet = 0;
    double *dxi  = new double[num_particles * Nd];
    std::fill_n(dxi, num_particles * Nd, 1);
	cout<<"dxi "<<dxi[0]<<endl;
	double xi_current = -1.0;
	double xi[1];
    /* Initial conditions */
    double* mom_init    = new double[num_particles * Nbeads * Nd];
    double* x_init      = new double[num_particles * Nbeads * Nd];
	double* x_meas_init = new double[num_particles * Nbeads * Nd];
    double* mom_minus_init    = new double[num_particles * Nbeads * Nd];
		
    double *centroid = new double[num_particles * Nd];


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
    int num_trajectories = 20;




    double corr[7000] = {0};
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
			constraint = 1;
            e_vel_verlet += velocity_verlet_integrate(thermostat_freq, mom_init, x_init, mom_out, x_out, sq_en, constraint, dxi, xi,  xi_current);
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
        double dum_energy;
//Long trajectory after thermlization
        constraint = 1;
        dum_energy = velocity_verlet_integrate(10000, mom_init, x_init, mom_out, x_out, sq_en, constraint, dxi, xi,  xi_current);

        double dum;
        //calculation step
        //x[0] = 0;
		double dev;
		int num_meas_traj = 50;
	    double pplus[num_meas_traj] = {0};
        double pminus[num_meas_traj] = {0};
		for(int i=0; i<num_particles; i++) {
			for(int j=0; j<Nbeads; j++) {
				for (int k = 0; k < Nd; k++) {
					*(x_meas_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_out + i * Nbeads * Nd  + j * Nd + k);
                    //cout<<"Position after verlet "<<*(x_meas_init + i * Nbeads * Nd  + j * Nd + k)<<endl;
				}
			}
		}
		
		//cout<<"I am here"<<endl;
		for(int meas_traj=1; meas_traj<=num_meas_traj;meas_traj++) {
			
			//pplus[0] = 0;
			//pminus[0] = 0;
        	for(int i=0; i<num_particles; i++) {
             
				dev =  sqrt(mass(i)/beta);
				for(int j=0; j<Nbeads; j++) {
                	for (int k = 0; k < Nd; k++) {
                    	std::normal_distribution<double> new_nd(0, dev);
						std::default_random_engine new_de(traj * meas_traj*(k + 1)*(j+10));//seed
						*(mom_init + i * Nbeads * Nd  + j * Nd + k) = new_nd(new_de);
						*(mom_minus_init + i * Nbeads * Nd  + j * Nd + k) =(-1)*(*(mom_init + i * Nbeads * Nd  + j * Nd + k));
						
                    	pplus[meas_traj] += *(mom_init + i * Nbeads * Nd  + j * Nd + k)/nbeads;
                    	pminus[meas_traj] += *(mom_minus_init + i * Nbeads * Nd  + j * Nd + k)/nbeads;

						*(x_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_meas_init + i * Nbeads * Nd  + j * Nd + k);
						
						
                	}
            	}	
        	}


        	for(int ch_t=1; ch_t<= 1200; ch_t++) {

	    		constraint = 0;
            	dum = velocity_verlet_integrate(1, mom_init, x_init, mom_out, x_out, sq_en, constraint, dxi, xi, xi_current);
               // cout<<"Energy plus "<<dum<<endl;			
				get_centroid(x_out,centroid);
				get_reaction_coordinates(centroid, xi_current, xi, dxi);
			
            	corr[ch_t] += heavside(xi[0]) * pplus[meas_traj];
                //cout<<"time "<<ch_t<<", plus "<< xi[0]<<endl;              	    	
            	for(int i=0; i<num_particles; i++) {
                	for(int j=0; j<Nbeads; j++) {
                    	for (int k = 0; k < Nd; k++) {
                            //cout<<"Position after "<<meas_traj<<" plus "<< *(x_out + i * Nbeads * Nd  + j * Nd + k)<<endl;
                        	*(x_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_out + i * Nbeads * Nd  + j * Nd + k);
                        	*(mom_init + i * Nbeads * Nd  + j * Nd + k) =  *(mom_out + i * Nbeads * Nd  + j * Nd + k);
                    	}
                	}
            	}
      		}
		
			for(int i=0; i<num_particles; i++) {
				for(int j=0; j<Nbeads; j++) {
					for (int k = 0; k < Nd; k++) {
						*(x_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_meas_init + i * Nbeads * Nd  + j * Nd + k);
					}
				}
			}

        	for(int ch_t=1; ch_t<= 1200; ch_t++) {

	    		constraint = 0;
            	dum = velocity_verlet_integrate(1, mom_minus_init, x_init, mom_out, x_out, sq_en, constraint, dxi, xi, xi_current);
                //cout<<"Energy minus "<<dum<<endl;
				get_centroid(x_out,centroid);
				get_reaction_coordinates(centroid, xi_current, xi, dxi);
			    //cout<<"time "<<ch_t<<", minus "<<corr[ch_t]<<endl; 
            	corr[ch_t] += heavside(xi[0]) * pminus[meas_traj];
			    //cout<<"time "<<ch_t<<", minus "<< xi[0]<<endl; 

                //cout<<"minus "<<xi[0]<<endl;
	    	    for(int i=0; i<num_particles; i++) {
                	for(int j=0; j<Nbeads; j++) {
                    	for (int k = 0; k < Nd; k++) {
                        	*(x_init + i * Nbeads * Nd  + j * Nd + k) =  *(x_out + i * Nbeads * Nd  + j * Nd + k);
                        	*(mom_minus_init + i * Nbeads * Nd  + j * Nd + k) =  *(mom_out + i * Nbeads * Nd  + j * Nd + k);
                    	}
                	}
            	}
      		}
	
		}
        cout<<"Traj number"<<traj<<"\t\t"<<corr[1]<<endl;
    }
   	for (int ch_t=0; ch_t <=1200;){

        vel<< ch_t*delt << "\t\t"<< corr[ch_t]/corr[1]<<endl;
        ch_t += 1;
    }

    //system("gnuplot plot 'vel_verlet_fluc.dat'");
	delete [] centroid;
	delete [] x_meas_init;
	delete [] mom_minus_init;
	delete [] dxi;
	delete [] x_init;
	delete [] x_out;
	delete [] mom_init;
	delete [] mom_out;
}

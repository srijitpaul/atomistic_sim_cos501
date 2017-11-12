#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
using namespace std;

double force_i(int i_particle, double x_i, double* x_j);
double etot(double* mom_j);
double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double *mom_out, double *x_out, double *sq_en);
double mass(int index);

double delt         =   0.01;
auto num_particles   =   1;
auto Nd              =   1;
auto Nbeads          =   1;
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
double force_i( int i_particle, double x_i, double* x_j){

    double force = 0;
    for(int j = 0; j<num_particles; j++) {

        if(j != i_particle)
            force += - mass(j) * (x_j[j] - x_i);
        else{}
    }
    if(num_particles==1)
        force = - mass(0) * x_i - 0.3 * pow(x_i,2) - 0.04 * pow(x_i,3);
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

double velocity_verlet_integrate(int npart, double* mom_init, double* x_init, double* mom_out, double* x_out, double *sq_en ) {

    double avg_en=0;
    double x_dum, mom_dum;
    double *x   = new double[num_particles * Nbeads * Nd];
    double *mom = new double[num_particles * Nbeads * Nd];



    /*Initial conditions loop*/
    for(int i=0; i<num_particles; i++) {
        for(int j=0; j<Nbeads; j++) {
            for (int k = 0; k < Nd; k++) {
                *(x + i * Nbeads * Nd  + j * Nd + k) = *(x_init + i * Nbeads * Nd + j * Nd + k);
                *(mom + i * Nbeads * Nd  + j * Nd + k) = *(mom_init + i * Nbeads * Nd + j * Nd + k);

            }
        }
    }
    avg_en += etot(mom);


    for(int i=0; i<=(npart   -   1); i++) {
            for (int p_index = 0; p_index < num_particles; p_index++) {
                for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
                    for (int dim = 0; dim < Nd; dim++) {

                        x_dum = *(x + p_index * Nbeads * Nd + bead_index * Nd + dim) +
                                    (*(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) / mass(p_index)) * delt +
                                    0.5 * (force_i(p_index, *(x + p_index * Nbeads * Nd + bead_index * Nd + dim), x)/ mass(p_index))
                                    * pow(delt, 2);

                        mom_dum = *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) +
                                        0.5 * (force_i(p_index, *(x + p_index * Nbeads * Nd + dim), x) +
                                                force_i(p_index, x_dum, x) )/ (mass(p_index)) * delt;

                        *(x + p_index * Nbeads * Nd + bead_index * Nd + dim) = x_dum;
                        *(mom + p_index * Nbeads * Nd + bead_index * Nd + dim) = mom_dum;

                        avg_en += etot(mom);
                        sq_en[0] += pow(etot(mom),2);
                    }

                }
            }
    }
    for (int p_index = 0; p_index < num_particles; p_index++) {
        for (int bead_index = 0; bead_index < Nbeads; bead_index++) {
            for (int dim = 0; dim < Nd; dim++) {
                *(x_out + p_index * Nbeads * Nd + bead_index * Nd + dim) = *(x + p_index * Nbeads * Nd +
                                                                             bead_index * Nd + dim);
                *(mom_out + p_index * Nbeads * Nd + bead_index * Nd + dim) = *(mom + p_index * Nbeads * Nd +
                                                                               bead_index * Nd + dim);
            }
        }
    }
    //cout<<*(x_out +1)<<endl;
    delete [] x;
    delete [] mom;
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

    //std::ofstream vel_verlet("vel_verlet_etot.dat", ios::out);
    std::ofstream vel("vel_verlet_fluc.dat", ios::out);
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
    int num_trajectories = 2000;
    double beta = 1;
    double x[60000] = {0};
    double mom[60000] = {0};
    double corr[30000] = {0};
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

                sigma[i] = sqrt(mass(i)/beta);
                std::normal_distribution<double> nd(0, sigma[i]);

                for(int j=0; j<Nbeads; j++) {
                    for (int k = 0; k < Nd; k++) {
                        std::default_random_engine de(traj*steps);//seed
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
        double dum;
        //calculation step
        x[0] = 0;
        mom[0] = 0;
        for(int i=0; i<num_particles; i++) {
            for(int j=0; j<Nbeads; j++) {
                for (int k = 0; k < Nd; k++) {
                    x[0] += *(x_out + i * Nbeads * Nd  + j * Nd + k);
                    mom[0] += *(mom_out + i * Nbeads * Nd  + j * Nd + k) ;
                    if(j==1) cout<<x[0]<<endl;
                }
            }
        }
        corr[0] += x[0] * x[0];

        for(int ch_t=1; ch_t<= 3000; ch_t++) {

            x[ch_t] = 0;
            mom[ch_t] = 0;
            dum = velocity_verlet_integrate(1, mom_init, x_init, mom_out, x_out, sq_en);
            for(int i=0; i<num_particles; i++) {
                for(int j=0; j<Nbeads; j++) {
                    for (int k = 0; k < Nd; k++) {
                         x[ch_t] += *(x_out + i * Nbeads * Nd  + j * Nd + k);
                         mom[ch_t] += *(mom_out + i * Nbeads * Nd  + j * Nd + k);
                    }
                }
            }
            corr[ch_t] += x[0] * x[ch_t];
            //cout<<x[0]<<endl;
            x_init = x_out;
            mom_init = mom_out;

        }
        cout<<"Traj number"<<traj<<"\t\t"<<dum<<endl;

        //vel << pow(1/beta,2)<<"\t\t"<<delt * total_steps * thermostat_freq << "\t\t"
        //          << sq_e_vel_verlet/(total_steps * thermostat_freq) - pow(e_vel_verlet / (total_steps * thermostat_freq),2) <<endl;
    }
    for (int ch_t=0; ch_t <=3000; ch_t++){

        vel<< ch_t*delt << "\t\t"<< corr[ch_t]/num_trajectories<<endl;
    }

    //system("gnuplot plot 'vel_verlet_fluc.dat'");
}
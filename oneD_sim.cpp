/////////////////////////////////////////////////
/////////   1-D simulation  /////////////////////
/////////////////////////////////////////////////

#include "psopt.h"
#include <math.h>

#define M2 13.6//8 // kg
#define M1S 14 //8.5 //kg
#define MF_INIT 0.6 //kg
#define G 9.81 // m/s2
#define THRUST_D 20 // N
#define THRUST_U 15
#define FLOW_RATE 0.018 // kg/s

#define TOTAL_NODE 100

using namespace PSOPT;

////// Define the end point cost function ///////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    //return tf;
    return 0;

}

///////// Define Lagrange cost function /////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    adouble xt_error = 0;
    // Approximate square wave function
    adouble x_ref, y, f;
    f = 1/(20);
    //y = (1/atan(1/0.001))*atan(sin(2*pi*time*f)/0.001);
    //x_ref = (y + smooth_fabs(y, 0.01))/2;
    //x_ref = sin(2*pi*time*f);
    if (time < 8)
    {
        x_ref = 1;
    }
    
    if (time >= 8)
    {
        x_ref = 0.1;
    }

    xt_error = states[0] - x_ref;

    
    //return 0.5*input*input;
    return 2*pow(xt_error,2);
    //return  0.0;
}

/////// Define system dynamics ///////
void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls, 
         adouble* patameters, adouble& time, adouble* xad, int iphase, Workspace* workspace)
         {
            adouble xdot, vdot, mdot, fddot, fudot, fdout, fuout;

            adouble x = states[0];
            adouble v = states[1];
            adouble m = states[2];
            //adouble fd = states[3];  //TRY
            //adouble fu = states[4];  //TRY

            adouble u1 = controls[0];
            adouble u2 = controls[1];

            ////////////TRY ///////////////////
            /*fddot = -27.7081*fd + 1.9349*u1;
            fdout = 215.9683*fd;

            fudot = -27.7081*fu + 1.9349*u2;
            fuout = 215.9683*fu;*/
            /////////////TRY ////////////////
            xdot = v;
            vdot = (M2-M1S-m)*G/(M2+M1S+m) + THRUST_D*u1/*fdout*//(M2+M1S+m)
                 - THRUST_U*u2/*fuout*//(M2+M1S+m);
            mdot = -FLOW_RATE * (u1+u2);//*smooth_fabs(u, 0.01);

            derivatives[0] = xdot;
            derivatives[1] = vdot;
            derivatives[2] = mdot;
            //derivatives[3] = fddot;// TRY
            //erivatives[4] = fudot; //TRY

            //path[0] = (u-1)*(u+1)*u ; //TODO: There is a problem here.          
         }

////////// Define the event function //////////
void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    adouble x0 = initial_states[0];
    adouble v0 = initial_states[1];
    adouble m0 = initial_states[2];
    //adouble fd0 = initial_states[3];
    //adouble fu0 = initial_states[4];
    //adouble xf = final_states[0];
    //adouble vf = final_states[1];

    e[0] = x0;
    e[1] = v0;
    e[2] = m0;
    //e[3] = xf;
    //e[4] = vf;
    //e[5] = fd0;
    //e[6] = fu0;

    
    
    
}

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    

}

////////// Define main routine //////////
int main(void)
{

    // declare key structures //
    Alg algorithm;
    Sol solution;
    Prob problem;
    Workspace ws(problem, algorithm, solution);  // TODO A PROBLEM HERE

    printf("CONFIG DONE!");

    // Rigister problem name //
    problem.name = "1-D Simulation";
    problem.outfilename = "1-d_sim.txt";

    // Define problem level constants & do level 1 setup //
    problem.nphases = 1;
    problem.nlinkages = 1;

    psopt_level1_setup(problem);

    // Define phase related information & do level 2 setup //
    problem.phases(1).nstates = 3;  //TRY
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents = 3;//5;
    problem.phases(1).npath = 0;// 1;
    problem.phases(1).nodes << TOTAL_NODE;

    psopt_level2_setup(problem, algorithm);

    ////////// Problem bounds information //////////
    problem.phases(1).bounds.lower.states << 0, -2, 0;//, 0, 0;  //TRY
    problem.phases(1).bounds.upper.states << 3, 2, MF_INIT;//, 15, 15; //TRY

    problem.phases(1).bounds.lower.controls << 0.0, 0.0;
    problem.phases(1).bounds.upper.controls << 1.0, 1.0;

    //problem.phases(1).bounds.lower.path(0) = 0.0;
    //problem.phases(1).bounds.upper.path(0) = 0.0;

    problem.phases(1).bounds.lower.events << 0, 0, MF_INIT;//, 0.7, 0;//, 0, 0;//TRY
    problem.phases(1).bounds.upper.events << 0, 0, MF_INIT;//, 0.7, 0;//, 0 ,0;// TRY

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = 10;
    problem.phases(1).bounds.upper.EndTime = 10;//2.5;

    ////////// Rigister problem functions //////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae;
    problem.events = &events;
    problem.linkages = &linkages;

    ////////// Initial guess //////////
    int nnodes= problem.phases(1).nodes(0);

    MatrixXd x0(3 /*5*/,nnodes); //TRY

    x0.row(0) = linspace(0.1,1.0, nnodes);
    x0.row(1) = linspace(0.1,1.0, nnodes);
    x0.row(2) = linspace(0.1,1.0, nnodes);
    //x0.row(3) = linspace(0.1,1.0, nnodes); //TRY
    //x0.row(4) = linspace(0.1,1.0, nnodes); //TRY

    problem.phases(1).guess.controls       = zeros(2,nnodes);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 2.0, nnodes);

    ////////// Algorithm options //////////
    algorithm.nlp_method = "IPOPT";
    algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic";
    algorithm.nlp_iter_max = 300;
    algorithm.nlp_tolerance = 1e-4;

    algorithm.collocation_method = "Legendre";

    

    
    ////////// Call PSOPT to solve the problem //////////
    psopt(solution, problem, algorithm);
    

    if (solution.error_flag) exit(0);
      

    ////////// Extract variables from solution structure //////////
    MatrixXd xStar 		= solution.get_states_in_phase(1);
    MatrixXd uStar 		= solution.get_controls_in_phase(1);
    MatrixXd t 		= solution.get_time_in_phase(1);
    MatrixXd H           = solution.get_dual_hamiltonian_in_phase(1);
    MatrixXd lambda      = solution.get_dual_costates_in_phase(1);

    /// Test section ///
    adouble interp_control;
    MatrixXd interp_state;
    //interp_state.resize(3,1);
    MatrixXd xStar_height = xStar.row(0);
    MatrixXd interp_time;
    interp_time.resize(1,10);
    for(int i = 0; i<10; i++)
    {
        interp_time(0,i) = 0.1*i;
    }
    //adouble interp_time = 0.5;
    //adouble& ref_time = interp_time;
     //spline_interpolation (&interp_control, ref_time, t, uStar, nnodes);
    printf("start to interpolate~~ \n");
    lagrange_interpolation(interp_state, interp_time, t, xStar_height);
    printf("finish! \n");
    //printf(" time %fs, height = %fm, velocity = %f m/s \n", interp_time(0,0), interp_state(0,0), interp_state(1,0));
    //lagrange_interpolation_ad(&interp_control, ref_time, )
    //printf("control @ time %fs: %f\n",interp_time.value(), interp_control.value());



    ////////// Save Solution data //////////
    Save(xStar,"x.dat");
    Save(uStar,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"p.dat");

    



    ////////// Plot the result //////////
    plot(t,xStar.block<3,TOTAL_NODE>(0,0),problem.name + ": states", "time (s)", "states", "x v m "); //TRY
    //plot(t,xStar.block<2,TOTAL_NODE>(3,0),problem.name + ": states", "time (s)", "states", "fd fu ");

    plot(interp_time, interp_state, problem.name + ": interp states", "time (s)", "interp states", "x"); // FOR INTERPOLATE TEST

    plot(t,uStar,problem.name + ": control", "time (s)", "control", "fd fu");

    plot(t,H,problem.name + ": Hamiltonian", "time (s)", "H", "H");

    plot(t,lambda,problem.name + ": costates", "time (s)", "lambda", "lambda_1 lambda_2 lambda_3");

    plot(t,xStar,problem.name + ": states", "time (s)", "states", "x y v",
                                  "pdf", "brac1_states.pdf");

    plot(t,uStar,problem.name + ": control", "time (s)", "control", "u",
                              "pdf", "brac1_control.pdf");

    plot(t,H,problem.name + ": Hamiltonian", "time (s)", "H", "H",
                                  "pdf", "brac1_hamiltonian.pdf");

    plot(t,lambda,problem.name + ": costates", "time (s)", "lambda", "lambda_1 lambda_2 lambda_3",
                                  "pdf", "brac1_costates.pdf");




}
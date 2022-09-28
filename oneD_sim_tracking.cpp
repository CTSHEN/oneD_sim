/////////////////////////////////////////////////
/////////   1-D simulation  /////////////////////
/////////////////////////////////////////////////

#include "psopt.h"
#include <math.h>

#define M2 9 // kg
#define M1S 8.5 //kg
#define MF_INIT 0.5 //kg
#define G 9.81 // m/s2
#define THRUST 15 // N
#define FLOW_RATE 0.018 // kg/s

#define TOTAL_NODE 50

#define DIS_X  1.5
#define DIS_V  0

using namespace PSOPT;

////// Define the end point cost function ///////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    // TODO: add end point error cost
    //adouble xf_error = final_states[0] - DIS_X;
    //adouble vf_error = final_states[1] - DIS_V;

    //return 1*(xf_error*xf_error + vf_error*vf_error);
    return 0;


}

///////// Define Lagrange cost function /////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    // TODO: Add control square cost
    //adouble input = controls[0];
    adouble xt_error = 0;
    // Approximate square wave function
    adouble x_ref, y, f;
    f = 1/(2*pi);
    //y = (1/atan(1/0.001))*atan(sin(2*pi*time*f)/0.001);
    //x_ref = (y + smooth_fabs(y, 0.01))/2;
    x_ref = sin(2*pi*time*f);

    xt_error = states[0] - x_ref;

    
    //return 0.5*input*input;
    return pow(xt_error,2);
}

/////// Define system dynamics ///////
void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls, 
         adouble* patameters, adouble& time, adouble* xad, int iphase, Workspace* workspace)
         {
            adouble xdot, vdot, mdot;

            adouble x = states[0];
            adouble v = states[1];
            adouble m = states[2];

            adouble u = controls[0];

            xdot = v;
            vdot = (M2-M1S-m)*G/(M2+M1S+m) + THRUST*u/(M2+M1S+m);
            mdot = -FLOW_RATE*smooth_fabs(u, 0.01);

            derivatives[0] = xdot;
            derivatives[1] = vdot;
            derivatives[2] = mdot;                   
         }

////////// Define the event function //////////
void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    adouble x0 = initial_states[0];
    adouble v0 = initial_states[1];
    adouble m0 = initial_states[2];
    //adouble xf = final_states[0];
    //adouble vf = final_states[1];

    e[0] = x0;
    e[1] = v0;
    e[2] = m0;
    //e[3] = xf;
    //e[4] = vf;
    //TODO cancel final state constraint?;
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

    // Rigister problem name //
    problem.name = "1-D Simulation";
    problem.outfilename = "1-d_sim.txt";

    // Define problem level constants & do level 1 setup //
    problem.nphases = 1;
    problem.nlinkages = 0;

    psopt_level1_setup(problem);

    // Define phase related information & do level 2 setup //
    problem.phases(1).nstates = 3;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents = 3;
    problem.phases(1).npath = 0;
    problem.phases(1).nodes << TOTAL_NODE;

    psopt_level2_setup(problem, algorithm);

    ////////// Problem bounds information //////////
    problem.phases(1).bounds.lower.states << 0, 0, 0;
    problem.phases(1).bounds.upper.states << 3, 2, MF_INIT;

    problem.phases(1).bounds.lower.controls << -1.0;
    problem.phases(1).bounds.upper.controls << 1.0;

    problem.phases(1).bounds.lower.events << 0, 0, MF_INIT;//, 1.4, 0;
    problem.phases(1).bounds.upper.events << 0, 0, MF_INIT;//, 1.6, 0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = 7.0;
    problem.phases(1).bounds.upper.EndTime = 8;

    ////////// Rigister problem functions //////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae;
    problem.events = &events;
    problem.linkages = &linkages;

    ////////// Initial guess //////////
    int nnodes= problem.phases(1).nodes(0);

    MatrixXd x0(3,nnodes);

    x0.row(0) = linspace(0.0,1.0, nnodes);
    x0.row(1) = linspace(0.0,1.0, nnodes);
    x0.row(2) = linspace(0.0,1.0, nnodes);

    problem.phases(1).guess.controls       = zeros(1,nnodes);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 2.0, nnodes);

    ////////// Algorithm options //////////
    algorithm.nlp_method = "IPOPT";
    algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic";
    algorithm.nlp_iter_max = 1000;
    algorithm.nlp_tolerance = 1e-6;

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
    adouble interp_time = 0.1;
    adouble& ref_time = interp_time;
    spline_interpolation (&interp_control, ref_time, t, uStar, nnodes);
    //get_interpolated_control(&interp_control, 0,1,ref_time, xad, workspace);
    //lagrange_interpolation_ad(&interp_control, ref_time, )
    printf("control @ time %fs: %f\n",interp_time.value(), interp_control.value());


    ////////// Save Solution data //////////
    Save(xStar,"x.dat");
    Save(uStar,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"p.dat");


    ////////// Plot the result //////////
    plot(t,xStar,problem.name + ": states", "time (s)", "states", "x v m");

    plot(t,uStar,problem.name + ": control", "time (s)", "control", "u");

    plot(t,H,problem.name + ": Hamiltonian", "time (s)", "H", "H");

    plot(t,lambda,problem.name + ": costates", "time (s)", "lambda", "lambda_1 lambda_2 lambda_3");

    plot(t,xStar,problem.name + ": states", "time (s)", "states", "x v m",
                                  "pdf", "brac1_states.pdf");

    plot(t,uStar,problem.name + ": control", "time (s)", "control", "u",
                              "pdf", "brac1_control.pdf");

    plot(t,H,problem.name + ": Hamiltonian", "time (s)", "H", "H",
                                  "pdf", "brac1_hamiltonian.pdf");

    plot(t,lambda,problem.name + ": costates", "time (s)", "lambda", "lambda_1 lambda_2 lambda_3",
                                  "pdf", "brac1_costates.pdf");




}
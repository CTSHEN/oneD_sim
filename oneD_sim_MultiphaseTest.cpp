/////////////////////////////////////////////////
/////   1-D simulation  with 2 phase test ///////
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

#define TOTAL_NODE 20
#define TOTAL_PHASE 8

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
    adouble x_ref;
    /*
    if (time < 9)
    {
        x_ref = 1;
    }
    
    if (time >= 9)
    {
        x_ref = 0.1;
    }*/
    x_ref = 1;

    xt_error = states[0] - x_ref;

    
    
    return 2*pow(xt_error,2);
    //return  0.0;
}

/////// Define system dynamics ///////
void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls, 
         adouble* patameters, adouble& time, adouble* xad, int iphase, Workspace* workspace)
         {
            adouble xdot, vdot, mdot, fddot, fudot, fdout, fuout;

            adouble x =     states[0];
            adouble v =     states[1];
            adouble m =     states[2];
            
            
            adouble u1 = controls[0];
            adouble u2 = controls[1];

            // dynamic of phase 1: only ThrusterD
            if (iphase == 1 || iphase==3 || iphase==5 || iphase==7)
            {
                xdot = v;
                vdot = (M2-M1S-m)*G/(M2+M1S+m) + THRUST_D*u1/(M2+M1S+m) -
                 THRUST_U*u2/(M2+M1S+m);
                
                mdot = -FLOW_RATE*(u1+u2);
            }

            // dynamic of phase 2 : only ThrusterU
            if(iphase == 2 || iphase==4 || iphase==6 || iphase==8)
            {
                xdot = v;
                vdot = (M2-M1S-m)*G/(M2+M1S+m) + THRUST_D*u1/(M2+M1S+m)-
                 THRUST_U*u1/(M2+M1S+m);
                mdot = -FLOW_RATE*(u1+u2);

            }
            

            derivatives[0] = xdot;
            derivatives[1] = vdot;
            derivatives[2] = mdot;
            
                
         }

////////// Define the event function //////////
void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    if (iphase == 1)
    {
        adouble x0 =    initial_states[0];
        adouble v0 =    initial_states[1];
        adouble m0 =    initial_states[2];
        

        e[0] = x0;
        e[1] = v0;
        e[2] = m0;
        
    }

    adouble xf2[5];
    if(iphase ==2)
    {
        get_final_states(   xf2, xad, 2, workspace);
        adouble vf2 = xf2[1];

        e[0] = vf2;

    }
    
    if (iphase == 8)
    {
        adouble xf = final_states[0];  //final state for phase 2
        adouble vf = final_states[1];  //final state for phase 2

        //e[0] = xf;
        //e[0] = vf;

    }


    
    
}

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    // Number of linkages:
    // Boundary of phase 1 to 8: 3 state continuity + 1 time continuity
    // Total: 28 linkage constraints

    adouble xf1[3], xi2[3], tf1, ti2; // P1 & P2
    adouble xf2[3], xi3[3], tf2, ti3; // P2 & P3
    adouble xf3[3], xi4[3], tf3, ti4; // P3 & P4
    adouble xf4[3], xi5[3], tf4, ti5; // P4 & P5
    adouble xf5[3], xi6[3], tf5, ti6; // P5 & P6
    adouble xf6[3], xi7[3], tf6, ti7; // P6 & P7
    adouble xf7[3], xi8[3], tf7, ti8; // P7 & P8

    // Linking phase 1 and 2
    get_final_states(   xf1, xad, 1, workspace);
    get_initial_states( xi2, xad, 2, workspace);
    tf1 = get_final_time(    xad, 1, workspace);
    ti2 = get_initial_time(  xad, 2, workspace);

    linkages[0] = xf1[0] - xi2[0];
    linkages[1] = xf1[1] - xi2[1];
    linkages[2] = xf1[2] - xi2[2];
    linkages[3] = tf1 - ti2;

    // Linking phase 2 and 3
    get_final_states(   xf2, xad, 2, workspace);
    get_initial_states( xi3, xad, 3, workspace);
    tf2 = get_final_time(    xad, 2, workspace);
    ti3 = get_initial_time(  xad, 3, workspace);

    linkages[4] = xf2[0] - xi3[0];
    linkages[5] = xf2[1] - xi3[1];
    linkages[6] = xf2[2] - xi3[2];
    linkages[7] = tf2 - ti3;

    // Linking phase 3 and 4
    get_final_states(   xf3, xad, 3, workspace);
    get_initial_states( xi4, xad, 4, workspace);
    tf3 = get_final_time(    xad, 3, workspace);
    ti4 = get_initial_time(  xad, 4, workspace);

    linkages[8] = xf3[0] - xi4[0];
    linkages[9] = xf3[1] - xi4[1];
    linkages[10] = xf3[2] - xi4[2];
    linkages[11] = tf3 - ti4;

    // Linking phase 4 and 5
    get_final_states(   xf4, xad, 4, workspace);
    get_initial_states( xi5, xad, 5, workspace);
    tf4 = get_final_time(    xad, 4, workspace);
    ti5 = get_initial_time(  xad, 5, workspace);

    linkages[12] = xf4[0] - xi5[0];
    linkages[13] = xf4[1] - xi5[1];
    linkages[14] = xf4[2] - xi5[2];
    linkages[15] = tf4 - ti5;

    // Linking phase 5 and 6
    get_final_states(   xf5, xad, 5, workspace);
    get_initial_states( xi6, xad, 6, workspace);
    tf5 = get_final_time(    xad, 5, workspace);
    ti6 = get_initial_time(  xad, 6, workspace);

    linkages[16] = xf5[0] - xi6[0];
    linkages[17] = xf5[1] - xi6[1];
    linkages[18] = xf5[2] - xi6[2];
    linkages[19] = tf5 - ti6;

    // Linking phase 6 and 7
    get_final_states(   xf6, xad, 6, workspace);
    get_initial_states( xi7, xad, 7, workspace);
    tf6 = get_final_time(    xad, 6, workspace);
    ti7 = get_initial_time(  xad, 7, workspace);

    linkages[20] = xf6[0] - xi7[0];
    linkages[21] = xf6[1] - xi7[1];
    linkages[22] = xf6[2] - xi7[2];
    linkages[23] = tf6 - ti7;

    // Linking phase 7 and 8
    get_final_states(   xf7, xad, 7, workspace);
    get_initial_states( xi8, xad, 8, workspace);
    tf7 = get_final_time(    xad, 7, workspace);
    ti8 = get_initial_time(  xad, 8, workspace);

    linkages[24] = xf7[0] - xi8[0];
    linkages[25] = xf7[1] - xi8[1];
    linkages[26] = xf7[2] - xi8[2];
    linkages[27] = tf7 - ti8;
}

////////// Define main routine //////////
int main(void)
{
    int PhaseIdx = 0;
    
    // declare key structures //
    Alg algorithm;
    Sol solution;
    Prob problem;
    Workspace ws(problem, algorithm, solution);  // TODO A PROBLEM HERE

    printf("CONFIG DONE!");

    // Rigister problem name //
    problem.name = "1-D Simulation Multuphase";
    problem.outfilename = "1-d_sim_Multuphase.txt";

    // Define problem level constants & do level 1 setup //
    problem.nphases = TOTAL_PHASE;
    problem.nlinkages = 28;

    psopt_level1_setup(problem);

    // Define phase related information & do level 2 setup //

    problem.phases(1).nstates = 3;
    problem.phases(1).ncontrols = 2;  
    problem.phases(1).nevents = 3;
    problem.phases(1).npath = 0;
    problem.phases(1).nodes << TOTAL_NODE;

    problem.phases(2).nstates = 3;
    problem.phases(2).ncontrols = 2;  
    problem.phases(2).nevents = 1;
    problem.phases(2).npath = 0;
    problem.phases(2).nodes << TOTAL_NODE;

    for(PhaseIdx=3;PhaseIdx<8; PhaseIdx++ )
    {
        problem.phases(PhaseIdx).nstates = 3;
        problem.phases(PhaseIdx).ncontrols = 2;  
        problem.phases(PhaseIdx).nevents = 0;
        problem.phases(PhaseIdx).npath = 0;
        problem.phases(PhaseIdx).nodes << TOTAL_NODE;
    }

    problem.phases(8).nstates = 3;
    problem.phases(8).ncontrols = 2;  
    problem.phases(8).nevents = 0;//2;
    problem.phases(8).npath = 0;
    problem.phases(8).nodes << TOTAL_NODE;

    psopt_level2_setup(problem, algorithm);

    ////////// Problem bounds information //////////
    // Phase 1
    problem.phases(1).bounds.lower.states << 0, -2, 0;
    problem.phases(1).bounds.upper.states << 3, 2, MF_INIT;

    problem.phases(1).bounds.lower.controls << 1.0, 0.0;
    problem.phases(1).bounds.upper.controls << 1.0, 0.0;

    problem.phases(1).bounds.lower.events << 0, 0, MF_INIT;
    problem.phases(1).bounds.upper.events << 0, 0, MF_INIT;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = 0;
    problem.phases(1).bounds.upper.EndTime = 8;

    for(PhaseIdx=3; PhaseIdx<=7; PhaseIdx = PhaseIdx+2)
    {
        problem.phases(PhaseIdx).bounds.lower.states << 0, -2, 0;
        problem.phases(PhaseIdx).bounds.upper.states << 3, 2, MF_INIT;
        problem.phases(PhaseIdx).bounds.lower.controls << 1.0, 0.0;
        problem.phases(PhaseIdx).bounds.upper.controls << 1.0, 0.0;
        problem.phases(PhaseIdx).bounds.lower.StartTime = 0.0;
        problem.phases(PhaseIdx).bounds.upper.StartTime = 8.0;
        problem.phases(PhaseIdx).bounds.lower.EndTime = 0;
        problem.phases(PhaseIdx).bounds.upper.EndTime = 8;
    }

    for(PhaseIdx=2; PhaseIdx<=6; PhaseIdx = PhaseIdx+2)
    {
        problem.phases(PhaseIdx).bounds.lower.states << 0, -2, 0;
        problem.phases(PhaseIdx).bounds.upper.states << 3, 2, MF_INIT;
        problem.phases(PhaseIdx).bounds.lower.controls << 0.0, 1.0;
        problem.phases(PhaseIdx).bounds.upper.controls << 0.0, 1.0;
        problem.phases(PhaseIdx).bounds.lower.StartTime = 0.0;
        problem.phases(PhaseIdx).bounds.upper.StartTime = 8.0;
        problem.phases(PhaseIdx).bounds.lower.EndTime = 0;
        problem.phases(PhaseIdx).bounds.upper.EndTime = 8;
    }

    problem.phases(2).bounds.lower.events << 0;
    problem.phases(2).bounds.upper.events << 0;

    // Phase 8
    problem.phases(8).bounds.lower.states << 0, -2, 0;
    problem.phases(8).bounds.upper.states << 3, 2, MF_INIT;

    problem.phases(8).bounds.lower.controls << 0.0, 1.0;
    problem.phases(8).bounds.upper.controls << 0.0, 1.0;

    //problem.phases(8).bounds.lower.events << -0.05;
    //problem.phases(8).bounds.upper.events << 0.05;

    problem.phases(8).bounds.lower.StartTime = 0.0;
    problem.phases(8).bounds.upper.StartTime = 8.0;

    problem.phases(8).bounds.lower.EndTime = 8;
    problem.phases(8).bounds.upper.EndTime = 8;

    ////////// Rigister problem functions //////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae;
    problem.events = &events;
    problem.linkages = &linkages;

    ////////// Initial guess //////////
    int nnodes= problem.phases(1).nodes(0);

    MatrixXd x0(3, nnodes); //TRY

    x0.row(0) = linspace(0.1,1.0, nnodes);
    x0.row(1) = linspace(0.1,1.0, nnodes);
    x0.row(2) = linspace(0.6,0.1, nnodes);
    

    //problem.phases(1).guess.controls       = zeros(2,nnodes);
    for(PhaseIdx=1; PhaseIdx<=8; PhaseIdx++)
    {
        problem.phases(PhaseIdx).guess.states         = x0;
        problem.phases(PhaseIdx).guess.time           = linspace(0.0, 8.0, nnodes);
    }
   
    

    ////////// Algorithm options //////////
    algorithm.nlp_method = "IPOPT";
    algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic";
    algorithm.nlp_iter_max = 300;
    algorithm.nlp_tolerance = 1e-3;

    algorithm.collocation_method = "Legendre";

    

    
    ////////// Call PSOPT to solve the problem //////////
    psopt(solution, problem, algorithm);
    

    if (solution.error_flag) exit(0);
      

    ////////// Extract variables from solution structure //////////
    MatrixXd xPhase1 		= solution.get_states_in_phase(1);
    MatrixXd uPhase1        = solution.get_controls_in_phase(1);
    MatrixXd tPhase1 		= solution.get_time_in_phase(1);
    //MatrixXd HPhase1           = solution.get_dual_hamiltonian_in_phase(1);
    //MatrixXd lambdaPhase1      = solution.get_dual_costates_in_phase(1);

    MatrixXd xPhase2 		= solution.get_states_in_phase(2);
    MatrixXd uPhase2        = solution.get_controls_in_phase(2);
    MatrixXd tPhase2 		= solution.get_time_in_phase(2);
    //MatrixXd HPhase2           = solution.get_dual_hamiltonian_in_phase(2);
    //MatrixXd lambdaPhase2      = solution.get_dual_costates_in_phase(2);
    MatrixXd xPhase3 		= solution.get_states_in_phase(3);
    MatrixXd uPhase3        = solution.get_controls_in_phase(3);
    MatrixXd tPhase3 		= solution.get_time_in_phase(3);

    MatrixXd xPhase4 		= solution.get_states_in_phase(4);
    MatrixXd uPhase4        = solution.get_controls_in_phase(4);
    MatrixXd tPhase4 		= solution.get_time_in_phase(4);

    MatrixXd xPhase5 		= solution.get_states_in_phase(5);
    MatrixXd uPhase5        = solution.get_controls_in_phase(5);
    MatrixXd tPhase5 		= solution.get_time_in_phase(5);

    MatrixXd xPhase6 		= solution.get_states_in_phase(6);
    MatrixXd uPhase6        = solution.get_controls_in_phase(6);
    MatrixXd tPhase6 		= solution.get_time_in_phase(6);

    MatrixXd xPhase7 		= solution.get_states_in_phase(7);
    MatrixXd uPhase7        = solution.get_controls_in_phase(7);
    MatrixXd tPhase7 		= solution.get_time_in_phase(7);

    MatrixXd xPhase8 		= solution.get_states_in_phase(8);
    MatrixXd uPhase8        = solution.get_controls_in_phase(8);
    MatrixXd tPhase8 		= solution.get_time_in_phase(8);

    MatrixXd xStar(3, length(tPhase1)+length(tPhase2)+length(tPhase3)
        +length(tPhase4)+length(tPhase5)+length(tPhase6)+length(tPhase7)+length(tPhase8));
    xStar << xPhase1, xPhase2, xPhase3, xPhase4, xPhase5, xPhase6, xPhase7, xPhase8;
    MatrixXd uStar(2, length(tPhase1)+length(tPhase2)+length(tPhase3)
        +length(tPhase4)+length(tPhase5)+length(tPhase6)+length(tPhase7)+length(tPhase8));
    uStar << uPhase1, uPhase2, uPhase3, uPhase4, uPhase5, uPhase6, uPhase7, uPhase8;
    MatrixXd tStar(1, length(tPhase1)+length(tPhase2)+length(tPhase3)
        +length(tPhase4)+length(tPhase5)+length(tPhase6)+length(tPhase7)+length(tPhase8));
    tStar << tPhase1, tPhase2, tPhase3, tPhase4, tPhase5, tPhase6, tPhase7, tPhase8;


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
    
    printf("start to interpolate~~ \n");
    lagrange_interpolation(interp_state, interp_time, tStar, xStar_height);
    printf("finish! \n");
    


    ////////// Save Solution data //////////
    //Save(xStar,"x.dat");
    //Save(uStar,"u.dat");
    //Save(t,"t.dat");
    //Save(lambda,"p.dat");

    



    ////////// Plot the result //////////
    plot(tStar,xStar.block<3,TOTAL_NODE*TOTAL_PHASE>(0,0),problem.name + ": states", "time (s)", "states", "x v m ");
    //plot(interp_time, interp_state, problem.name + ": interp states", "time (s)", "interp states", "x"); // FOR INTERPOLATE TEST

    multiplot(tStar,uStar,problem.name + ": control", "time (s)", "control", "ThrusterDown ThrusterUp", 2, 1);

    //plot(t,H,problem.name + ": Hamiltonian", "time (s)", "H", "H");

    //plot(t,lambda,problem.name + ": costates", "time (s)", "lambda", "lambda_1 lambda_2 lambda_3");

    



}

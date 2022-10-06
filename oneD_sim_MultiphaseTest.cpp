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
#define TOTAL_PHASE 14
#define FINAL_TIME 10

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
    
    /*if (time < 5)
    {
        x_ref = 1;
    }
    
    if (time >= 5)
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
            //if (iphase == 1 || iphase==3 || iphase==5 || iphase==7 || iphase==9 || iphase==11)
            if(iphase%2 != 0)
            {
                xdot = v;
                vdot = (M2-M1S-m)*G/(M2+M1S+m) + THRUST_D*u1/(M2+M1S+m) -
                 THRUST_U*u2/(M2+M1S+m);
                
                mdot = -FLOW_RATE*(u1+u2);
            }

            // dynamic of phase 2 : only ThrusterU
            //if(iphase == 2 || iphase==4 || iphase==6 || iphase==8 || iphase==10 || iphase==12)
            if(iphase%2 == 0)
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

        //e[0] = vf2;

    }
    
    if (iphase == TOTAL_PHASE)
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
    // Boundary of phase 1 to 12: 3 state continuity + 1 time continuity
    // Total: 44 linkage constraints

    int index = 0;

    auto_link(linkages, &index, xad, 1, 2,  workspace);
    auto_link(linkages, &index, xad, 2, 3,  workspace);
    auto_link(linkages, &index, xad, 3, 4,  workspace);
    auto_link(linkages, &index, xad, 4, 5,  workspace);
    auto_link(linkages, &index, xad, 5, 6,  workspace);
    auto_link(linkages, &index, xad, 6, 7,  workspace);
    auto_link(linkages, &index, xad, 7, 8,  workspace);
    auto_link(linkages, &index, xad, 8, 9,  workspace);
    auto_link(linkages, &index, xad, 9, 10, workspace);
    auto_link(linkages, &index, xad, 10,11, workspace);
    auto_link(linkages, &index, xad, 11, 12, workspace);
    auto_link(linkages, &index, xad, 12,13, workspace);
    auto_link(linkages, &index, xad, 13, 14, workspace);
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
    problem.nlinkages = 4*(TOTAL_PHASE-1);

    psopt_level1_setup(problem);

    // Define phase related information & do level 2 setup //

    problem.phases(1).nstates = 3;
    problem.phases(1).ncontrols = 2;  
    problem.phases(1).nevents = 3;
    problem.phases(1).npath = 0;
    problem.phases(1).nodes << TOTAL_NODE;

    problem.phases(2).nstates = 3;
    problem.phases(2).ncontrols = 2;  
    problem.phases(2).nevents = 0;
    problem.phases(2).npath = 0;
    problem.phases(2).nodes << TOTAL_NODE;

    for(PhaseIdx=3;PhaseIdx<TOTAL_PHASE; PhaseIdx++ )
    {
        problem.phases(PhaseIdx).nstates = 3;
        problem.phases(PhaseIdx).ncontrols = 2;  
        problem.phases(PhaseIdx).nevents = 0;
        problem.phases(PhaseIdx).npath = 0;
        problem.phases(PhaseIdx).nodes << TOTAL_NODE;
    }

    problem.phases(TOTAL_PHASE).nstates = 3;
    problem.phases(TOTAL_PHASE).ncontrols = 2;  
    problem.phases(TOTAL_PHASE).nevents = 0;//2;
    problem.phases(TOTAL_PHASE).npath = 0;
    problem.phases(TOTAL_PHASE).nodes << TOTAL_NODE;

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
    problem.phases(1).bounds.upper.EndTime = FINAL_TIME;

    for(PhaseIdx=3; PhaseIdx<=TOTAL_PHASE; PhaseIdx = PhaseIdx+2)
    {
        problem.phases(PhaseIdx).bounds.lower.states << 0, -2, 0;
        problem.phases(PhaseIdx).bounds.upper.states << 3, 2, MF_INIT;
        problem.phases(PhaseIdx).bounds.lower.controls << 1.0, 0.0;
        problem.phases(PhaseIdx).bounds.upper.controls << 1.0, 0.0;
        problem.phases(PhaseIdx).bounds.lower.StartTime = 0.0;
        problem.phases(PhaseIdx).bounds.upper.StartTime = FINAL_TIME;
        problem.phases(PhaseIdx).bounds.lower.EndTime = 0;
        problem.phases(PhaseIdx).bounds.upper.EndTime = FINAL_TIME;
    }

    for(PhaseIdx=2; PhaseIdx<TOTAL_PHASE; PhaseIdx = PhaseIdx+2)
    {
        problem.phases(PhaseIdx).bounds.lower.states << 0, -2, 0;
        problem.phases(PhaseIdx).bounds.upper.states << 3, 2, MF_INIT;
        problem.phases(PhaseIdx).bounds.lower.controls << 0.0, 1.0;
        problem.phases(PhaseIdx).bounds.upper.controls << 0.0, 1.0;
        problem.phases(PhaseIdx).bounds.lower.StartTime = 0.0;
        problem.phases(PhaseIdx).bounds.upper.StartTime = FINAL_TIME;
        problem.phases(PhaseIdx).bounds.lower.EndTime = 0;
        problem.phases(PhaseIdx).bounds.upper.EndTime = FINAL_TIME;
    }

    //problem.phases(2).bounds.lower.events << 0;
    //problem.phases(2).bounds.upper.events << 0;

    // Last phase
    problem.phases(TOTAL_PHASE).bounds.lower.states << 0, -2, 0;
    problem.phases(TOTAL_PHASE).bounds.upper.states << 3, 2, MF_INIT;
    problem.phases(TOTAL_PHASE).bounds.lower.controls << 0.0, 1.0;
    problem.phases(TOTAL_PHASE).bounds.upper.controls << 0.0, 1.0;
    problem.phases(TOTAL_PHASE).bounds.lower.StartTime = 0.0;
    problem.phases(TOTAL_PHASE).bounds.upper.StartTime = FINAL_TIME;
    problem.phases(TOTAL_PHASE).bounds.lower.EndTime = FINAL_TIME;
    problem.phases(TOTAL_PHASE).bounds.upper.EndTime = FINAL_TIME;

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
    for(PhaseIdx=1; PhaseIdx<=TOTAL_PHASE; PhaseIdx++)
    {
        problem.phases(PhaseIdx).guess.states         = x0;
        problem.phases(PhaseIdx).guess.time           = linspace(0.0, FINAL_TIME, nnodes);
    }
   
    

    ////////// Algorithm options //////////
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 300;
    algorithm.nlp_tolerance               = 1e-3;
    algorithm.mesh_refinement             = "automatic";
    algorithm.collocation_method          = "Legendre";

    

    
    ////////// Call PSOPT to solve the problem //////////
    psopt(solution, problem, algorithm);
    

    if (solution.error_flag) exit(0);
      
    MatrixXd xStar(3, TOTAL_PHASE*TOTAL_NODE);
    MatrixXd uStar(2, TOTAL_PHASE*TOTAL_NODE);
    MatrixXd tStar(1, TOTAL_PHASE*TOTAL_NODE);

    for(PhaseIdx=1; PhaseIdx<=TOTAL_PHASE; PhaseIdx++ )
    {
        xStar.block<3,TOTAL_NODE>(0,TOTAL_NODE*(PhaseIdx-1)) = solution.get_states_in_phase(PhaseIdx);
        uStar.block<2,TOTAL_NODE>(0,TOTAL_NODE*(PhaseIdx-1)) = solution.get_controls_in_phase(PhaseIdx);
        tStar.block<1,TOTAL_NODE>(0,TOTAL_NODE*(PhaseIdx-1)) = solution.get_time_in_phase(PhaseIdx);
    }
       


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

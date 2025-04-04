
#ifndef _CSPGEN3
#define _CSPGEN3

#include <cmath>

using namespace std;

class cspGen3CostModel {
public:

    cspGen3CostModel();             // Class initializer. Memory is shared and modified between optimizations. 
    ~cspGen3CostModel() = default;  // Default destructor. Memory is out-of-scope after the optimization is completed. 
    void designRoutine();           // Accepts a power block to design a CSP plant and calculate LCOE. 

    // derived parameters    
    double W_dot_therm;     // [MWt]    power cycle thermal input
    double W_dot_field;     // [MWt]    power delivered to the receiver from the solar field
    double W_dot_rec;       // [MWt]    power required at and absorbed by the receiver
    double W_dot_less;      // [MWe]    net power, less parasitics
    double W_dot_losses;    // [MWt]    receiver estimated losses
    double W_elec_annual;   // [MWe-h]  annual electricity produced

    struct cycle { // power block data structure

        // power block parameters
        double T_phx_i;     // [K]   pHX particle inlet temperature
        double T_phx_o;     // [K]   pHX particle outlet temperature
        double T_HTR_i;     // [K]   HTR sCO2 inlet temperature (cold side)
        double T_HTR_o;     // [K]   HTR sCO2 outlet temperature (cold side)
        double T_LTR_i;     // [K]   LTR sCO2 inlet temperature (cold side)
        double T_LTR_o;     // [K]   LTR sCO2 outlet temperature (cold side)
        double T_trb_i;     // [K]   turbine inlet temperature
        double eta_gen;     // [-]   electrical generator efficiency
        double efficiency;  // [-]   cycle efficiency
        double W_dot_net;   // [MWt] cycle design power output (W_dot_t - W_dot_mc - W_dot_rc)
        double W_dot_gen;   // [MWe] cycle net power generation (W_dot_net * eta_gen)
        double phx_height;  // [m]   primary heat exchanger height
        double P_max;       // [MPa] power cycle high pressure
        double P_min;       // [MPa] power cycle low pressure

        cycle() {
            T_phx_i = T_phx_o = eta_gen = efficiency = W_dot_net = T_HTR_i = T_HTR_o = 
                T_LTR_i = T_LTR_o = W_dot_gen = phx_height = T_trb_i = P_max = P_min = 0.0;
        };
    } s_cycle;

    struct costs { // capital costs data structure

        // CSP capital costs
        double solar_tower;                 // [$] CSP Gen3 Solar Tower capital cost
        double solar_field;                 // [$] CSP Gen3 Solar Field capital cost
        double falling_particle_receiver;   // [$] CSP Gen3 Falling Particle Receiver capital cost
        double particles;                   // [$] bulk cost of particles
        double particle_losses;             // [$] incurred cost due to particle loss / attrition
        double particle_storage;            // [$] particle storage bins and insulation capital cost
        double particle_lifts;              // [$] particle transportation capital cost
        double land;                        // [$] bulk cost of land required for power plant

        // Cycle capital costs
        double HTR_capital;            // [$] high temperature recuperator capital cost
        double LTR_capital;            // [$] low temperature recuperator capital cost
        double PHX_capital;            // [$] primary heat exchanger capital cost
        double air_cooler_capital;     // [$] air cooler capital cost
        double compressor_capital;     // [$] primary compressor capital cost
        double recompressor_capital;   // [$] recompressor capital cost
        double turbine_capital;        // [$] turbine capital cost
        double piping_inventory_etc;   // [$] piping, inventory control, etc. 
        double balance_of_plant;       // [$] transformers, inverters, controls, etc.

        // Total capital, maintenance, and LCOE
        double cycle_capital;               // [$]       Power block capital costs
        double plant_capital;               // [$]       CSP equipment capital costs
        double total_capital;               // [$]       total expected capital cost of plant
        double annual_maintenance;          // [$/year]  expected O&M annual costs
        double total_adjusted_cost;         // [$]       total cost of capital, construction, and contingencies 
        double levelized_cost_of_energy;    // [$/MWe-h] CSP Gen3 levelized cost of energy

        costs() {
            solar_tower = solar_field = falling_particle_receiver = land = balance_of_plant = 
                particles = particle_storage = particle_lifts = particle_losses = 0;

            HTR_capital = LTR_capital = PHX_capital = air_cooler_capital =
                compressor_capital = recompressor_capital = turbine_capital = 0;

            annual_maintenance = total_adjusted_cost = cycle_capital = piping_inventory_etc = 
                plant_capital = total_capital = levelized_cost_of_energy = 0;
        };
    } s_costs;

    struct financing { // financing data structure

        // assumed lifetime, maintenance, and financing rates
        double discount_rate = 0.070;    // [-] cost of financing         (Albrecht, 2019)
        double construction = 0.060;     // [-] construction costs        (Albrecht, 2019)
        double contingency = 0.100;      // [-] unexpected costs          (Albrecht, 2019)
        double indirect = 0.130;         // [-] indirect cost of capital  (Albrecht, 2019)
        double inflation = 0.025;        // [-] average rate of inflation (DOE, SunShot Vision Study)
        double lifetime = 30.00;         // [years] total design lifetime (DOE, SunShot Vision Study)
        double maintenance = 40000;      // [$/MWe-year] O&M costs        (Albrecht, 2019)
        double balance_of_plant = 167E3; // [$/MWe] balance of plant rate (NREL, SAM) 

        // financing derived parameters
        double real_discount_rate;       // [-] calculated using discount rate and inflation
        double capital_recovery_factor;  // [-] calculated using real discount rate and assumed lifetime

        financing() {
            real_discount_rate = ((1 + discount_rate) / (1 + inflation)) - 1; // [-]
            capital_recovery_factor =
                real_discount_rate * pow(1 + real_discount_rate, lifetime)
                / (pow(1 + real_discount_rate, lifetime) - 1);
        };
    } s_financing;

    struct parasitics { // net power parasitics data structure

        // net power parasitics
        double field;   // [MWe] parasitics from heliostat tracking
        double lifts;   // [MWe] parasitics from particle lifts
        double cooler;  // [MWe] parasitics from power block air cooler

        parasitics() {
            field = lifts = cooler = 0;
        };
    } s_parasitics;

    struct storage { // thermal energy storage data structure

        // thermal energy storage capacity
        double hours_of_capacity;   // [hours]  hours of energy storage
        double capacity_factor;     // [-]      actual energy delivered / nameplate
        double dT;                  // [K]      storage losses

        struct warm { // warm storage

            // warm particle storage dimensions and performance
            double height;      // [m]
            double radius;      // [m]
            double volume;      // [m2]
            double efficiency;  // [-]
            double Tm;          // [K] mean temperature
            double Ti;          // [K] inlet temperature
            double To;          // [K] outlet temperature

            warm() {
                height = radius = volume = efficiency = 0;
                Tm = Ti = To = 0; 
            }
        } s_warm;

        struct cold { // cold storage

            // cold particle storage dimensions and performance
            double height;      // [m]
            double radius;      // [m]
            double volume;      // [m2]
            double efficiency;  // [-]
            double Tm;          // [K] mean temperature
            double Ti;          // [K] inlet temperature
            double To;          // [K] outlet temperature

            cold() {
                height = radius = volume = efficiency = 0;
                Tm = Ti = To = 0;
            }
        } s_cold;

        storage() {
            hours_of_capacity = capacity_factor = dT = 0;
        };
    } s_storage;

    struct particles { // thermal carrier data structure

        // parameters derived from particle
        double cost_per_kg;     // [$/kg]   particle bulk cost
        double angle_of_repose; // [rad]    particle repose angle
        double bulk_density;    // [kg/m3]  particle bulk density

        // parameters derived from cycle
        double m_particles;     // [kg]     bulk particle mass
        double m_dot_phx;       // [kg/s]   mass flow through phx
        double m_dot_rec;       // [kg/s]   mass flow through receiver

        // assumed non-storage particle factor
        double non_storage = 0.05; // [-] 

        particles() {
            cost_per_kg = angle_of_repose = bulk_density
                = m_particles = m_dot_phx = m_dot_rec = 0;
        };
    } s_particles;

    struct field { // solar field data structure

        // ratio of receiver thermal input to power cycle thermal input
        double solar_multiple;  // [-] solar multiple determines hours of energy storage

        // assumed solar field parameters
        double heliostat_cost_per_area = 75; // [$/m2]  cost per unit area of heliostat           (US DOE, 2012)
        double prep_cost_per_area = 10.0;    // [$/m2]  cost per unit area of land preparation    (US DOE, 2012)
        double land_cost_per_area = 2.47;    // [$/m2]  cost per unit area of land                (Mehos, 2016)
        double tracking_power = 0.0055;      // [kW/kW] heliostat tracking power gross fraction   (NREL SAM, 2024)
        double heliostat_x = 12.0;           // [m]     heliostat height
        double heliostat_y = 12.0;           // [m]     heliostat width

        // derived solar field parameters
        double area_heliostats; // [m2] total heliostat surface area required
        double area_total_land; // [m2] total land required

        field() {
            area_heliostats = area_total_land = solar_multiple = 0;
        };
    } s_field;

    struct tower { // solar tower data structure

        // solar tower dimensions
        double height; // [m]   tower height
        double radius; // [m]   tower radius
        double Vb;     // [m/s] wind velocity at base of tower
        double Vt;     // [m/s] wind velocity at top of tower

        tower() {
            height = radius = Vt = Vb = 0;
        };
    } s_tower;

    struct receiver { // falling particle receiver data structure

        // falling particle receiver dimensions and performance
        double height;               // [m]   flat plate receiver height 
        double width;                // [m]   flat plate receiver width
        double wind_direction;       // [deg] wind direction (North=0, North-facing receiver)
        double area_aperature;       // [m2]  aperature area = receiver area
        double aspect_ratio;         // [-]   receiver height / receiver width
        double particle_loss_factor; // [-]   particle loss from open air receiver
        double efficiency_modifier;  // [-]   receiver efficiency modifier
        double efficiency;           // [-]   receiver efficiency
        double Tm;                   // [K]   particle mean temperature in receiver 
        double Ti;                   // [K]   particle temperature at receiver inlet
        double To;                   // [K]   particle temperature at receiver outlet 

        receiver() {
            height = width = area_aperature = aspect_ratio = wind_direction =
                particle_loss_factor = efficiency = Tm =
                Ti = To = efficiency_modifier = 0.0;
        };
    } s_receiver;

    struct lifts { // particle transportation data structure

        // particle transportation / lift parameters
        double height;      // [m] total combined lift height
        double efficiency;  // [-] lift electrical efficiency
        double dT;          // [K] losses in lift

        lifts() {
            height = efficiency = dT = 0;
        };
    } s_lifts;

    struct piping { // power block piping data structure

        // piping design parameters and requirements
        double lifetime;                // [hours]  total cumulative hours of operation before failure
        double stress_creep_rupture;    // [MPa]    stress to creep rupture for the lifetime and temperature
        double thickness_ratio;         // [-]      ratio of the pipe wall thickness to outer diameter (th / ro)
        double baseline;                // [$/m]    cost per length of piping at 700C
        double cost_per_length;         // [$/m]    calculated cost per length of piping using the M-R-M parameterization
        double normalized_cost;         // [-]      normalized cost per length of piping
        double dP;                      // [MPa]    pressure difference across the wall of the pipe
        double ri;                      // [m]      pipe internal radius to meet mass flow / velocity assumptions
        double ro;                      // [m]      pipe external radius to meet stress-to-creep-rupture requirement
        double th;                      // [m]      pipe wall thickness
        double Ac;                      // [m2]     cross-sectional area of the pipe
        double T_max;                   // [C]      piping maximum temperature
        double f_fix;                   // [-]      piping cost factor, fixed component
        double f_var;                   // [-]      piping cost factor, variable component
        double cost_factor;             // [-]      piping cost factor, piping cost = this * power block cost

        struct alloy { // pipe alloy structure

            // physical and economic properties of selected alloy
            double specific_cost;   // [$/kg]   cost per kilogram of the selected alloy
            double density;         // [kg/m3]  density of the selected alloy
            string name;            // [-]      name of the selected alloy

            alloy() {
                specific_cost = density = 0.0;
                name = "NA"; 
            }
        } s_alloy;
        
        piping() {
            lifetime = dP = stress_creep_rupture = thickness_ratio = baseline = T_max = 
                cost_per_length = normalized_cost = ri = ro = th = Ac = cost_factor = 0.0;
            f_fix = 0.05; // [-] selected to fit data from White, et al. 2017. 
            f_var = 0.07; // [-] selected to fit data from White, et al. 2017.  
        }
    } s_piping;

private:

    enum { // Equipment types for CEPCI cost correction
        BASE = 0,   // Default (average CEPCI)
        PIPE = 1,   // Piping, Valves, Fittings, etc.
        HTEX,       // Heat Exchangers and Tanks
        TURB,       // Pumps and Turbomachinery
        LABR,       // Labor and Supervision
        LAND,       // Land and Buildings
        LIFT,       // Lifts / Hoists / Process Machinery
    };

    void temperatureCostScaling();  // Scales the cost of the LTR, HTR, and turbine based on the respective operating temperature 
    void particleTemperatures();    // Calculates particle temperatures at the inlet / outlet of the receiver and warm / cold storage
    void receiverLosses();          // Calculates estimated receiver thermal losses
    void sizeEquipment();           // Sizes CSP Gen3 equipment (solar tower, etc). 
    void pipingFactor();            // Calculates piping % share of power block costs.
    double costLand();              // Calculates cost of the total land required.
    double costTower();             // Calculates cost of the solar tower. 
    double costField();             // Calculates cost of the solar field / heliostats. 
    double costLifts();             // Calculates cost of the particle lifts. 
    double costStorage();           // Calculates cost of the thermal energy storage. 
    double costReceiver();          // Calculates cost of the falling particle receiver. 
    double costParticles();         // Calculates cost of bulk particles required. 
    double costPiping(double T);    // Calculates the piping cost per length for a given temperature 
    double costParticleLosses();    // Calculates incurred cost of particle loss / attrition. 

    double CEPCI(                   // Using the US BLS CPI, returns a correction factor for inflation to Jan 2025.
        int year, int type
    );   

};


#endif // _CSPGEN3


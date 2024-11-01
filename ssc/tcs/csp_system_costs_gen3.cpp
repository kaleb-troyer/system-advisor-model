
#include <cmath>
#include "csp_system_costs_gen3.h"

const double pi = 3.14159265358979323846;

cspGen3CostModel::cspGen3CostModel() {

    // initializing unstructured parameters
    W_dot_therm = W_dot_rec = W_dot_less = W_elec_annual = W_dot_field = W_dot_losses = 0;

    // initializing structs
    s_parasitics = parasitics();
    s_financing = financing();
    s_particles = particles();
    s_receiver = receiver();
    s_storage = storage();
    s_cycle = cycle();
    s_costs = costs();
    s_field = field();
    s_tower = tower();
    s_lifts = lifts();

};

void cspGen3CostModel::designRoutine() {
    /*
    Kaleb Troyer, 2024.

    CSP Gen3 design routine for calculating the levelized
    cost of energy for a given power block.

    Design Routine:
    1) Extract power cycle parameters.
    2) Estimate thermal energy storage hours and capacity factor.
    3) Estimate power required at the falling particle receiver.
    4) Retreive particulate properties and characteristics.
    5) Size all CSP Gen3 equipment.
    6) Calculate all estimated capital costs.
    7) Calculate estimated net annual energy generation.
    8) Calculate the levelized cost of energy.
    */

    // particulate properties / characteristics 
    s_particles.angle_of_repose = 0.559;
    s_particles.bulk_density = 1625;
    s_particles.cost_per_kg = 0.185;

    // solar multiple, thermal energy storage
    s_field.solar_multiple = 2.5;       // [-] (Albrecht, 2019) 
    s_storage.hours_of_capacity = 14;   // [h] (Albrecht, 2019) 
    s_storage.capacity_factor = 0.71;   // [-] (Albrecht, 2019) 

    // estimate efficiencies here
    s_storage.s_warm.efficiency = 0.98;     // assumed
    s_storage.s_cold.efficiency = 0.98;     // assumed
    s_receiver.efficiency = 0.9;            // assumed
    s_lifts.efficiency = 0.8;               // assumed

    // calculating power required at the receiver, according to cycle efficiency
    W_dot_therm = s_cycle.W_dot_net / s_cycle.efficiency; // [MWt] power required at power block
    W_dot_rec = s_field.solar_multiple * W_dot_therm;     // [MWt] power absorbed by receiver

    temperatures();

    double tolerance = W_dot_rec * 1E-4;    // convergence criteria for receiver losses
    double iterations = 0;                  // iteration count for convergence
    double max_iters = 10;                  // maximum iterations allowed for convergence

    do {
        /*
        Iteratively solving for equipment sizes until thermal losses
        in the receiver has converged. 
        */

        W_dot_field = W_dot_rec + W_dot_losses;  // Update field power requirements
        sizeEquipment();   // Recalculate equipment sizing based on new field power
        receiverLosses();  // Recalculate receiver losses based on current setup

        iterations++;
    } while (abs(W_dot_field - (W_dot_rec + W_dot_losses)) > tolerance && iterations < max_iters);

    // getting CSP costs from cost models
    s_costs.solar_tower = costTower();                  // [$]
    s_costs.solar_field = costField();                  // [$]
    s_costs.falling_particle_receiver = costReceiver(); // [$]
    s_costs.particles = costParticles();                // [$]
    s_costs.particle_losses = costParticleLosses();     // [$]
    s_costs.particle_storage = costStorage();           // [$]
    s_costs.particle_lifts = costLifts();               // [$]
    s_costs.land = costLand();                          // [$]

    // calculating parasitics and annual electricity produced
    s_parasitics.field = s_field.tracking_power * W_dot_field;                                          // [MWe]
    s_parasitics.lifts = 1E-6 * s_particles.m_dot_rec * s_lifts.height * 9.80665 / s_lifts.efficiency;  // [MWe]
    W_dot_less = s_cycle.W_dot_net - (s_parasitics.field + s_parasitics.lifts + s_parasitics.cooler);   // [MWe]
    W_elec_annual = W_dot_less * s_storage.capacity_factor * (24 * 365);                                // [MWe-h / year]
    
    // calculating levelized cost of energy
    s_costs.balance_of_plant = s_financing.balance_of_plant * s_cycle.W_dot_net; // [$]
    s_costs.cycle_capital = s_costs.HTR_capital_cost + s_costs.LTR_capital_cost + s_costs.PHX_capital_cost + s_costs.air_cooler_capital_cost + s_costs.compressor_capital_cost + s_costs.recompressor_capital_cost + s_costs.turbine_capital_cost;
    s_costs.piping_inventory_etc = s_costs.cycle_capital * 0.2; 
    s_costs.cycle_capital += s_costs.piping_inventory_etc; 
    s_costs.plant_capital = s_costs.solar_tower + s_costs.solar_field + s_costs.falling_particle_receiver + s_costs.particles + s_costs.particle_losses + s_costs.particle_storage + s_costs.particle_lifts + s_costs.land + s_costs.balance_of_plant;
    s_costs.total_capital = s_costs.cycle_capital + s_costs.plant_capital; 
    s_costs.annual_maintenance = s_financing.maintenance * s_cycle.W_dot_net;
    s_costs.total_adjusted_cost = (1 + s_financing.construction) * (1 + s_financing.indirect) * (1 + s_financing.contingency) * s_costs.total_capital; 
    s_costs.levelized_cost_of_energy =
        ((s_financing.lifetime * s_costs.total_adjusted_cost * s_financing.capital_recovery_factor) + (s_financing.lifetime * s_costs.annual_maintenance))
        / (s_financing.lifetime * W_elec_annual);

};

void cspGen3CostModel::receiverLosses() {
    /*
    Falling-particle receiver thermal losses are estimated using the receiver
    temperature, aperature area, and power delivered to the receiver. Work done
    by Gonzales-Portillo (2020) is extended to cover solar field powers greater
    than 700MW and less than 375MW.

    Losses are then scaled to cover receiver temperatures above and below 700C.
    */

    const double sigma = 5.670374419E-8;              // Stefan-Boltzmann Constant 
    const double epsilon = 0.80;                      // receiver emissivity (Gonzales-Portillo, 2020) (assumed)
    const double heat_transfer_coefficient = 95;      // (Gonzales-Portillo, 2020) (assumed)
    const double T_baseline = 1073.15;                // (Gonzales-Portillo, 2020) (assumed)
    const double T_ambient = 308.15;                  // (Gonzales-Portillo, 2020) (assumed)

    // losses are scaled by receiver temperature
    double E1 = sigma * epsilon * (pow(s_receiver.T_avg, 4) - pow(T_ambient, 4)) + heat_transfer_coefficient * (s_receiver.T_avg - T_ambient);
    double E2 = sigma * epsilon * (pow(T_baseline, 4) - pow(T_ambient, 4)) + heat_transfer_coefficient * (T_baseline - T_ambient);
    double scale = E1 / E2; 

    // model developed from work done by Gonzales-Portillo, 2020. 
    double C1 = 2.771e-02; 
    double C2 = 5.245e-04; 
    double C3 = 6.403e-08; 
    double C4 = -5.735e-07; 
    W_dot_losses = (C1 + C2 * s_receiver.area_aperature
        + C3 * pow(s_receiver.area_aperature, 2)
        + C4 * s_receiver.area_aperature * W_dot_field) * W_dot_field * scale;

}; 

void cspGen3CostModel::sizeEquipment() {
    /*
    Linear regression models for CSP plant sizing were developed
    using SolarPILOT and scipy.optimize, by minimizing the total energy
    specific cost and fitting the regression to the top 10% of results.

    Input:
    -> Power delivered to Receiver [MW] (100 to 600)

    Output:
    -> Solar Tower [m]; RMSE = 12.980
    -> Solar Field [m^2]; RMSE = 32e3
    -> Receiver Width [m]; RMSE = 2.528
    */

    // Solar tower sizing
    s_tower.height = 4.821e+01 + (W_dot_field * 4.665e-01);
    s_tower.radius = 15;

    // Solar field sizing
    s_field.area_heliostats = -6.272e+04 + (W_dot_field * 2.174e+03);
    s_field.area_total_land = s_field.area_heliostats * 6 + 45000;

    // Receiver sizing
    s_receiver.height = 15;     // assumed in SolarPILOT study
    s_receiver.width = 9.113 + (3.274e-02 * W_dot_field);
    s_receiver.area_aperature = s_receiver.height * s_receiver.width;
    s_receiver.aspect_ratio = s_receiver.height / s_receiver.width;
    s_receiver.particle_loss_factor = 0.000001;     // assumed
    
    // Particle bulk calculations
    s_particles.m_particles = s_particles.m_dot_phx * s_storage.hours_of_capacity * 3600;
    s_particles.m_dot_rec = s_particles.m_dot_phx * s_field.solar_multiple; 

    // particle storage sizing
    double height_buffer = 2;       // [m] assumed
    s_storage.s_warm.radius = 12;   // [m] assumed
    s_storage.s_cold.radius = 12;   // [m] assumed
    s_storage.s_warm.volume = s_particles.m_particles / s_particles.bulk_density;
    s_storage.s_cold.volume = s_particles.m_particles / s_particles.bulk_density;
    s_storage.s_warm.height = height_buffer + (s_storage.s_warm.volume - (pi / 3) * pow(s_storage.s_warm.radius, 3) * tan(s_particles.angle_of_repose)) / (pi * pow(s_storage.s_warm.radius, 2));   // height from required volume of cylinder + cone, calculated with angle of repose
    s_storage.s_cold.height = height_buffer + (s_storage.s_cold.volume - (pi / 3) * pow(s_storage.s_cold.radius, 3) * tan(s_particles.angle_of_repose)) / (pi * pow(s_storage.s_cold.radius, 2));   // height from required volume of cylinder + cone, calculated with angle of repose

    // Determining if particle storage can fit inside tower to calculate total lift height
    if (s_cycle.phx_height == 0) { s_cycle.phx_height = s_storage.s_warm.height; }

    double tower_free_space = s_tower.height * 0.1;    // [m] space left intentionally free in tower (assumed)
    double height_required = s_storage.s_warm.height + s_storage.s_cold.height + s_cycle.phx_height + s_receiver.height;
    if (height_required < s_tower.height - tower_free_space) {
        s_lifts.height = height_required;
    }
    else {
        s_lifts.height = s_storage.s_warm.height + s_storage.s_cold.height + s_tower.height;
    }

};

void cspGen3CostModel::temperatures() {
    /*
    Estimating the receiver and TES inlet / outlet temperatures. These
    temperatures are derived from the PHX particle inlet temperature, which
    is a decision variable in the optimization process. 
    */

    double dT = 5; // TES temperature drop
    // thermal energy storage inlet / outlet temperatures
    s_storage.s_warm.To = s_cycle.T_phx_i; 
    s_storage.s_cold.Ti = s_cycle.T_phx_o; 
    s_storage.s_warm.Ti = s_storage.s_warm.To + dT; 
    s_storage.s_cold.To = s_storage.s_cold.Ti - dT; 
    s_storage.s_warm.T_avg = (s_storage.s_warm.Ti + s_storage.s_warm.To) / 2; 
    s_storage.s_cold.T_avg = (s_storage.s_cold.Ti + s_storage.s_cold.To) / 2; 

    // receiver particle inlet / outlet and body temperatures
    s_receiver.Ti = s_storage.s_cold.To; 
    s_receiver.To = s_storage.s_warm.Ti; 
    s_receiver.T_avg = (s_receiver.Ti + s_receiver.To + s_receiver.To + s_receiver.To) / 4; // taking a three-quarters average to estimate receiver temperature

}

double cspGen3CostModel::costTower() {
    /*
    R. Buck and S. Giuliano, "Impact of Solar Tower Design
    Parameters on sCO2-Based Solar Tower Plants," in 2nd
    European Supercritical CO2 Conference, Essen, 2018.
    */
    const double C1 = 157.440;  // [$/m^C2]
    const double C2 = 1.91744;  // [-]
    return C1 * pow(s_tower.height, C2);
};

double cspGen3CostModel::costField() {
    /*
    US DOE, "Sunshot vision study," Washington, DC, 2012.
    */
    return (s_field.heliostat_cost_per_area + s_field.prep_cost_per_area) * s_field.area_heliostats;
};

double cspGen3CostModel::costLifts() {
    /*
    K. D. Repole and S. M. Jeter, "Design and Analysis of a
    High Temperature Particulate Hoist for Proposed Particle
    Heating Concentrator Solar Power Systems," in ASME
    2016 10th International Conference on Energy
    Sustainability, Charlotte, 2016.
    */
    const double C1 = 58.37;    // [$-s/m-kg]
    return C1 * s_lifts.height * s_particles.m_dot_rec;
};

double cspGen3CostModel::costLand() {
    /*
    Mehos, M.; Turchi, C.; Jorgenson, J.; Denholm, P.; Ho, C.; Armijo, K. On the Path to SunShot: Advancing Concentrating Solar Power
    Technology, Performance, and Dispatchability; EERE Publication and Product Library: Golden, CO, USA, 2016.
    */
    return s_field.land_cost_per_area * s_field.area_total_land;
};

double cspGen3CostModel::costStorage() {
    /*
    Albrecht, K.J.; Bauer, M.L.; Ho, C.K. Parametric analysis of particle CSP system performance and cost to intrinsec particle
    properties and operating conditions. In Proceedings of the 13th International Conference on Energy and Sustainability, Bellevue,
    WA, USA, 14–17 July 2019.
    */
    const double C1 = 1230; // [$/m^2] 
    const double C2 = 0.37; // [$/m^2] 
    const double C3 = 600;  // [K] 
    const double C4 = 400;  // [K]
    double C_bin_warm = C1 + C2 * ((s_storage.s_warm.T_avg - C3) / C4); // [$/m2] cost of bin insulation
    double C_bin_cold = C1 + C2 * ((s_storage.s_cold.T_avg - C3) / C4); // [$/m2] cost of bin insulation
    double A_bin_warm = 2 * pi * s_storage.s_warm.radius * s_storage.s_warm.height + pi * s_storage.s_warm.radius * pow(pow(s_storage.s_warm.height, 2) + pow(s_storage.s_warm.radius, 2), 0.5);
    double A_bin_cold = 2 * pi * s_storage.s_cold.radius * s_storage.s_cold.height + pi * s_storage.s_cold.radius * pow(pow(s_storage.s_cold.height, 2) + pow(s_storage.s_cold.radius, 2), 0.5);
    return (C_bin_warm * A_bin_warm) + (C_bin_cold * A_bin_cold);
};

double cspGen3CostModel::costReceiver() {
    /*
    C. K. Ho, "A review of high-temperature particle receivers
    for concentrating solar power," Applied Thermal
    Engineering, vol. 109, pp. 958-969, 2016.
    */
    const double C1 = 37400; // [$/m^2]
    return C1 * s_receiver.height * s_receiver.width;
};

double cspGen3CostModel::costParticles() {
    /*
    Albrecht, K.J.; Bauer, M.L.; Ho, C.K. Parametric analysis of particle CSP system performance and cost to intrinsec particle
    properties and operating conditions. In Proceedings of the 13th International Conference on Energy and Sustainability, Bellevue,
    WA, USA, 14–17 July 2019.

    M. Mehos, C. Turchi, J. Vidal, M. Wagner, Z. Ma, C. Ho,
    W. Kolb, C. Andraka and A. Kruizenga, "Concentrating
    Solar Power Gen3 Demonstration Roadmap," National
    Renewable Energy Laboratory, Golden, 2017.
    */
    return (1 + s_particles.non_storage) * s_particles.cost_per_kg * s_particles.m_particles;
};

double cspGen3CostModel::costParticleLosses() {
    /*
    Albrecht, K.J.; Bauer, M.L.; Ho, C.K. Parametric analysis of particle CSP system performance and cost to intrinsec particle
    properties and operating conditions. In Proceedings of the 13th International Conference on Energy and Sustainability, Bellevue,
    WA, USA, 14–17 July 2019.
    */
    const double C1 = 365;  // [days / year]
    const double C2 = 3600; // [seconds / hour] 
    double annual_particle_losses = C1 * C2 * s_storage.hours_of_capacity * s_particles.m_dot_phx * s_receiver.particle_loss_factor;
    return s_financing.lifetime * annual_particle_losses * s_particles.cost_per_kg;
};



#include "cudaq/algorithms/evolve.h"
#include "cudaq/algorithms/integrator.h"
#include "cudaq/operators.h"
#include "export_csv_helper.h"
#include <cudaq.h>

template <size_t... Is>
void save_expectation_to_export_csv(const std::string& filename, 
                                    const std::vector<double>& tpoints, 
                                    auto& evolve_result,
                                    std::index_sequence<Is...>) 
{    
    auto all_exps = evolve_result.expectation_values.value();

    constexpr size_t num_spins = sizeof...(Is);
    const size_t num_time_steps = tpoints.size();
    
    std::vector<double> spin_columns[num_spins];
    
    for(size_t i = 0; i < num_spins; ++i) 
    {
        spin_columns[i].resize(num_time_steps);
    }

    for(size_t t = 0; t < num_time_steps; ++t) 
    {
        for(size_t i = 0; i < num_spins; ++i) 
        {
            spin_columns[i][t] = (double)all_exps[t][i];
        }
    }
    
    export_csv(filename, 
               {"time", tpoints}, 
               std::pair<std::string, const std::vector<double>&>{"sigma_z_" + std::to_string(Is), spin_columns[Is]}...);
}

int main() 
{
    auto start = std::chrono::steady_clock::now();
    
    const int64_t num_spins = 14;

    cudaq::dimension_map dimensions;

    for(int64_t i = 0; i < num_spins; ++i) 
    {
        dimensions[i] = 2;
    }

    std::string spin_state;

    for(int64_t i = num_spins - 1; i >= 0; --i) 
    {
        spin_state.push_back((i == 0) ? '1' : '0');
    }

    const int64_t initial_state_index = std::stoll(spin_state, nullptr, 2);
    const int64_t state_size = 1LL << num_spins;

    std::vector<std::complex<double>> psi0_data(state_size, {0.0, 0.0});
    psi0_data[initial_state_index] = {1.0, 0.0};

    auto psi0 = cudaq::state::from_data(psi0_data);

    const int64_t num_steps = 1001;
    std::vector<double> steps = cudaq::linspace(0.0, 100.0, num_steps);

    const double J = 1.0;

    std::vector<std::tuple<double, double, double>> positions = {{-4.89898, -4.89898, 9.79796}, {-5.49670, -2.66827, 11.4310}, 
                                                                 {-7.64733, -6.23311, 8.16497}, {-8.24504, -4.00240, 9.79796}, 
                                                                 {-8.84276, -1.77169, 11.4310}, {-6.23311, -7.64733, 8.16497}, 
                                                                 {-6.83083, -5.41662, 9.79796}, {-7.42855, -3.18591, 11.4310}, 
                                                                 {-8.02627, -0.95519, 13.0639}, {-5.41662, -6.83083, 9.79796}, 
                                                                 {-6.01433, -4.60012, 11.4310}, {-6.61205, -2.36941, 13.0639}, 
                                                                 {-4.00240, -8.24504, 9.79796}, {-4.60012, -6.01433, 11.4310}, 
                                                                 {-5.19784, -3.78362, 13.0639}, {-5.79555, -1.55291, 14.6969}, 
                                                                 {-3.18591, -7.42855, 11.4310}, {-3.78362, -5.19784, 13.0639}, 
                                                                 {-4.38134, -2.96713, 14.6969}, {-6.53197, -6.53197, 7.75672}, 
                                                                 {-7.12969, -4.30126, 9.38971}, {-7.72741, -2.07055, 11.0227}, 
                                                                 {-5.71548, -5.71548, 9.38971}, {-6.31319, -3.48477, 11.0227}, 
                                                                 {-6.91091, -1.25406, 12.6557}, {-4.30126, -7.12969, 9.38971}, 
                                                                 {-4.89898, -4.89898, 11.0227}, {-5.49670, -2.66827, 12.6557}, 
                                                                 {-3.48477, -6.31319, 11.0227}, {-4.08248, -4.08248, 12.6557}, 
                                                                 {-4.68020, -1.85177, 14.2887}, {-7.64733, -6.23311, 9.38971}, 
                                                                 {-8.24504, -4.00240, 11.0227}, {-8.84276, -1.77169, 12.6557}, 
                                                                 {-6.23311, -7.64733, 9.38971}, {-6.83083, -5.41662, 11.0227}, 
                                                                 {-7.42855, -3.18591, 12.6557}, {-5.41662, -6.83083, 11.0227}, 
                                                                 {-6.01433, -4.60012, 12.6557}, {-6.61205, -2.36941, 14.2887}, 
                                                                 {-4.00240, -8.24504, 11.0227}, {-4.60012, -6.01433, 12.6557}, 
                                                                 {-5.19784, -3.78362, 14.2887}};

    auto hamiltonian = cudaq::spin_op::empty();

    for(int64_t i = 0; i < num_spins; ++i) 
    {        
        for(int64_t j = i + 1; j < num_spins; ++j)
        {
            auto [xi, yi, zi] = positions[i];
            auto [xj, yj, zj] = positions[j];

            double dx = xj - xi;
            double dy = yj - yi;
            double dz = zj - zi;

            double dr3 = std::pow(dx * dx + dy * dy + dz * dz, 1.5);
            double dr = std::pow(dx * dx + dy * dy + dz * dz, 0.5);
            double cos_tetha = dz / dr;
            double Jzz = J * (1.0 - 3.0 * cos_tetha * cos_tetha) / dr3;
            
            hamiltonian += Jzz / 4.0 * cudaq::spin_op::z(i) * cudaq::spin_op::z(j);
            hamiltonian -= Jzz / 4.0 * (cudaq::spin_op::plus(i) * cudaq::spin_op::minus(j) + 
                                        cudaq::spin_op::plus(j) * cudaq::spin_op::minus(i));

        }
    }

    cudaq::schedule schedule(steps);

    cudaq::integrators::runge_kutta integrator(4, 0.1);

    std::vector<cudaq::spin_op> observables(num_spins);

    for(int64_t i = 0; i < num_spins; ++i)
    {
        observables[i] = cudaq::spin_op::z(i);
    }

    auto evolve_result = cudaq::evolve(hamiltonian, 
                                       dimensions, 
                                       schedule, 
                                       psi0, 
                                       integrator, 
                                       std::vector<cudaq::spin_op>{},
                                       observables,
                                       cudaq::IntermediateResultSave::ExpectationValue);

    save_expectation_to_export_csv("result.csv", 
                                   steps, 
                                   evolve_result, 
                                   std::make_index_sequence<num_spins>{});

    auto finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total spent time: " << elapsed.count() << " s." << std::endl;

    return 0;
}
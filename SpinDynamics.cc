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

std::vector<std::tuple<size_t, size_t, size_t>> SquareLattice(int64_t Nx, int64_t Ny)
{
    size_t Nbond = 2 * Nx * Ny;
    
    std::vector<std::tuple<size_t, size_t, size_t>> latt;
    latt.reserve(Nbond);

    auto mod = [](int64_t i, int64_t N)
    {
        return (i < 0) ? (i % N + N) % N : i % N;
    };
    
    for(int64_t i = 0; i < Nx; ++i)
    {
        for(int64_t j = 0; j < Ny; ++j)
        {
            latt.emplace_back(i + j * Nx, mod(i + 1, Nx) + j * Nx, 0);
            latt.emplace_back(i + j * Nx, i + mod(j + 1, Ny) * Nx, 1);
        } 
    }
    
    if(latt.size() != Nbond) std::cerr << "SquareLattice: Wrong number of bonds." << std::endl;
    
    return latt;
}

int main() 
{
    auto start = std::chrono::steady_clock::now();
    
    const int64_t num_spins = 25;

    cudaq::dimension_map dimensions;

    for(int64_t i = 0; i < num_spins; ++i) 
    {
        dimensions[i] = 2;
    }

    std::string spin_state;

    for(int64_t i = num_spins - 1; i >= 0; --i) 
    {
        spin_state.push_back('0');
    }

    const int64_t initial_state_index = std::stoll(spin_state, nullptr, 2);
    const int64_t state_size = 1LL << num_spins;

    std::vector<std::complex<double>> psi0_data(state_size, {0.0, 0.0});
    psi0_data[initial_state_index] = {1.0, 0.0};

    auto psi0 = cudaq::state::from_data(psi0_data);

    const int64_t num_steps = 201;
    std::vector<double> steps = cudaq::linspace(0.0, 1.0, num_steps);

    double hz = 3.04438;

    auto latt = SquareLattice(5L, 5L);

    auto hamiltonian = cudaq::spin_op::empty();

    for(int64_t i = 0; i < num_spins; ++i) 
    {        
        hamiltonian -= hz * cudaq::spin_op::z(i);
    }

    for(const auto& [i, j, b] : latt)
    {                                      
        hamiltonian -= cudaq::spin_op::x(i) * cudaq::spin_op::x(j);
    }

    cudaq::schedule schedule(steps);

    cudaq::integrators::magnus_expansion integrator(10, 0.005);

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
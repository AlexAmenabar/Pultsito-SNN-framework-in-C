#!/bin/bash



COMMANDS=(
"gcc -Wall -O3 -march=skylake-avx512 -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -g -o bin/main_simulation_rd src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DREORDER"
"gcc -Wall -O3 -march=skylake-avx512 -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation_rd src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512 -DREORDER"
)
#
#COMMANDS=(
#    "gcc -Wall -O2 -march=skylake-avx512 -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm"
#    "gcc -Wall -O2 -march=native -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_parallel_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c  -lm -DOPENMP -fopenmp"
#    "gcc -Wall -O2 -march=native -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_simulation_rd src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DREORDER"
#    "gcc -Wall -O2 -march=native -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_parallel_simulation_rd src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DOPENMP -fopenmp -DREORDER"
#    "gcc -Wall -O2 -march=native -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_parallel_samples_simulation_rd src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DOPENMP -fopenmp -DREORDER -DPAR_SAMPLES"
#    "gcc -Wall -O2 -march=native -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -g -o bin/main_parallel_samples_simulation_rd_debug src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DOPENMP -fopenmp -DREORDER -DPAR_SAMPLES"
#    "gcc -Wall -O2 -march=native -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_nested_parallel_simulation_rd src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DOPENMP -fopenmp -DREORDER -DPAR_SAMPLES -DNESTED"
#)

echo "== Compiling... =="

# ejecutar cada comando
for CMD in "${COMMANDS[@]}"; do
    echo "  Generating executable..."
    # ejecutar el comando
    eval $CMD
done

echo ""
echo "== Code compiled! =="
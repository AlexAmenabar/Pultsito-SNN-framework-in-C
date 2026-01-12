#!/bin/bash



COMMANDS=(
"gcc -O3 -ffast-math -fopenmp -Iinclude -Iinclude/pultsito -Ilib -g -o bin/main_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm"
"gcc -O3 -march=skylake-avx512 -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512"
"gcc -O3 -march=skylake-avx512 -ffast-math -funroll-loops -fprefetch-loop-arrays -flto -mtune=native -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation_opt src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512"
"gcc -O3 -ffast-math -fopenmp -funroll-loops -fprefetch-loop-arrays -flto -mtune=native -Iinclude -Iinclude/pultsito -Ilib -o bin/main_simulation_opt src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm"
)

#"gcc -O3 -march=native -ffast-math -funroll-loops -fprefetch-loop-arrays -flto -fprofile-generate -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512 -DREORDER"
#"gcc -O3 -march=native -ffast-math -funroll-loops -fprefetch-loop-arrays -flto --fprofile-use -fopenmp -Iinclude -Iinclude/pultsito -Ilib -o bin/main_avx_simulation src/pultsito/main.c src/pultsito/snn_library.c src/pultsito/load_data.c src/pultsito/helpers.c src/pultsito/neuron_models/lif_neuron.c src/pultsito/training_rules/stdp.c src/pultsito/simulations/simulations.c src/pultsito/metrics.c -lm -DAVX512 -DREORDER"

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